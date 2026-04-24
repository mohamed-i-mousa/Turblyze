/******************************************************************************
 * @file PvdTimeSeries.cpp
 * @brief Implementation of PVD collection file helpers
 *****************************************************************************/

// Own header
#include "PvdTimeSeries.h"

// STL includes
#include <fstream>
#include <iostream>
#include <vector>

// Header includes
#include "ErrorHandler.h"


namespace VTK
{

void writePVDTimeSeriesHeader
(
    const std::string& pvdFilename
)
{
    std::ofstream pvdFile(pvdFilename);

    if (!pvdFile.is_open())
    {
        FatalError("Failed to open PVD file: " + pvdFilename);
    }

    pvdFile << "<?xml version=\"1.0\"?>\n";
    pvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" "
            << "byte_order=\"LittleEndian\">\n";
    pvdFile << "  <Collection>\n";
    pvdFile << "  </Collection>\n";
    pvdFile << "</VTKFile>\n";

    pvdFile.close();

    std::cout
        << "PVD time series header written: "
        << pvdFilename << std::endl;
}

void appendPVDTimeStep
(
    const std::string& pvdFilename,
    const std::string& vtuFilename,
    Scalar timeValue
)
{
    // Read existing PVD file
    std::ifstream pvdFileIn(pvdFilename);
    if (!pvdFileIn.is_open())
    {
        FatalError
        (
            "Failed to open PVD file for reading: "
          + pvdFilename
        );
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(pvdFileIn, line))
    {
        lines.push_back(line);
    }
    pvdFileIn.close();

    // Find the </Collection> line and insert before it
    std::ofstream pvdFileOut(pvdFilename);
    if (!pvdFileOut.is_open())
    {
        FatalError
        (
            "Failed to open PVD file for writing: "
          + pvdFilename
        );
    }

    for (const auto& existingLine : lines)
    {
        if (existingLine.find("</Collection>") != std::string::npos)
        {
            // Insert the new timestep before closing collection
            pvdFileOut << "    <DataSet timestep=\"" << timeValue
                      << "\" file=\"" << vtuFilename << "\"/>\n";
        }
        pvdFileOut << existingLine << "\n";
    }

    pvdFileOut.close();

    std::cout
        << "Added timestep " << timeValue << " to PVD file"
        << std::endl;
}

} // namespace VTK
