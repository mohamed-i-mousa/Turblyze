/******************************************************************************
 * @file PvdTimeSeries.cpp
 * @brief Implementation of PVD collection file helpers
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "PvdTimeSeries.h"

// Standard library headers
#include <fstream>
#include <iostream>
#include <vector>

// Project headers
#include "ErrorHandler.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

void writePVDTimeSeriesHeader
(
    const FilePath& pvdFile
)
{
    std::ofstream pvdFileOutput(pvdFile);

    if (!pvdFileOutput.is_open())
    {
        FatalError("Failed to open PVD file: " + pvdFile);
    }

    pvdFileOutput
        << "<?xml version=\"1.0\"?>" << '\n'
        << "<VTKFile type=\"Collection\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << '\n'
        << "  <Collection>" << '\n'
        << "  </Collection>" << '\n'
        << "</VTKFile>" << '\n';

    pvdFileOutput.close();

    std::cout
        << "PVD time series header written: "
        << pvdFile << '\n';
}


void appendPVDTimeStep
(
    const FilePath& pvdFile,
    const FilePath& vtuFile,
    Scalar timeValue
)
{
    // Read existing PVD file
    std::ifstream pvdFileInput(pvdFile);
    if (!pvdFileInput.is_open())
    {
        FatalError
        (
            "Failed to open PVD file for reading: "
          + pvdFile
        );
    }

    std::vector<Message> lines;
    Message line;
    while (std::getline(pvdFileInput, line))
    {
        lines.push_back(line);
    }
    pvdFileInput.close();

    // Find the </Collection> line and insert before it
    std::ofstream pvdFileOutput(pvdFile);
    if (!pvdFileOutput.is_open())
    {
        FatalError
        (
            "Failed to open PVD file for writing: "
          + pvdFile
        );
    }

    bool inserted = false;

    for (const auto& existingLine : lines)
    {
        if (existingLine.find("</Collection>") != Message::npos)
        {
            // Insert the new timestep before closing collection
            pvdFileOutput
                << "    <DataSet timestep=\"" << timeValue
                << "\" file=\"" << vtuFile << "\"/>" << '\n';
            inserted = true;
        }
        pvdFileOutput << existingLine << '\n';
    }

    if (!inserted)
    {
        FatalError
        (
            "PVD file '" + pvdFile
          + "' has no </Collection> marker; cannot append timestep."
        );
    }

    pvdFileOutput.close();

    if (pvdFileOutput.fail())
    {
        FatalError("Failed to write PVD file: " + pvdFile);
    }

    std::cout
        << "Added timestep " << timeValue << " to PVD file" << '\n';
}

} // namespace VTK
