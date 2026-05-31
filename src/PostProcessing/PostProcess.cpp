/******************************************************************************
 * @file PostProcess.cpp
 * @brief After-solve reporting and result export
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "PostProcess.h"

// Standard library headers
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <string>
#include <string_view>

// Project headers
#include "CaseConfiguration.h"
#include "DerivedFields.h"
#include "ErrorHandler.h"
#include "Mesh.h"
#include "SIMPLE.h"
#include "VTK/VtkBoundaryWriter.h"
#include "VTK/VtkWriter.h"
#include "TurbulenceModel.h"

// *************************** namespace PostProcess **************************

namespace PostProcess
{

void reportStatistics
(
    const SIMPLE& solver,
    const CaseConfiguration& config
)
{
    std::cout
         << '\n' << "--- 5. Extracting Solution Fields ---" << '\n';

    const ScalarField& Ux = solver.Ux();
    const ScalarField& Uy = solver.Uy();
    const ScalarField& Uz = solver.Uz();
    const ScalarField& pressure = solver.pressure();

    if (config.debug)
    {
        std::cout
            << "Solution extracted." << '\n';
    }

    std::cout
         << '\n' << "--- 6. Post-Processing Results ---" << '\n';

    const ScalarField velocityMag = VTK::velocityMagnitude(Ux, Uy, Uz);

    if (Ux.empty())
    {
        Warning("Solution fields are empty. Skipping statistics.");
        return;
    }

    Scalar maximumVelocity = S(0.0);
    Scalar averageVelocity = S(0.0);
    Scalar maximumPressure = pressure[0];
    Scalar minimumPressure = pressure[0];

    for (size_t cellIdx = 0; cellIdx < Ux.size(); ++cellIdx)
    {
        const Scalar vmag = velocityMag[cellIdx];
        maximumVelocity = std::max(maximumVelocity, vmag);
        averageVelocity += vmag;

        maximumPressure = std::max(maximumPressure, pressure[cellIdx]);
        minimumPressure = std::min(minimumPressure, pressure[cellIdx]);
    }
    averageVelocity /= S(Ux.size());

    std::cout
        << "Flow Statistics:" << '\n';
    std::cout
        << "  Max velocity magnitude: " << maximumVelocity
        << " m/s" << '\n';
    std::cout
        << "  Average velocity magnitude: "
        << averageVelocity << " m/s" << '\n';
    std::cout
        << "  Pressure range: [" << minimumPressure
        << ", " << maximumPressure << "] Pa" << '\n';
}


void exportResults
(
    const SIMPLE& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const CaseConfiguration& config
)
{
    std::cout
         << '\n' << "--- 7. Exporting Results to VTK ---" << '\n';

    const ScalarField& Ux = solver.Ux();
    const ScalarField& Uy = solver.Uy();
    const ScalarField& Uz = solver.Uz();
    const ScalarField& pressure = solver.pressure();

    const ScalarField velocityMag = VTK::velocityMagnitude(Ux, Uy, Uz);

    std::map<std::string, const ScalarField*> scalarFieldsToVtk;

    scalarFieldsToVtk["pressure"] = &pressure;
    scalarFieldsToVtk["velocityMagnitude"] = &velocityMag;

    for (const auto& output : turbulence.cellDataOutputs())
    {
        if (output.second != nullptr)
        {
            scalarFieldsToVtk[std::string{output.first}] = output.second;
        }
    }

    std::map<std::string, std::array<const ScalarField*, 3>>
    vectorFieldsToVtk;

    vectorFieldsToVtk["velocity"] = {&Ux, &Uy, &Uz};

    std::string vtuFilename = config.vtkOutputFilename;
    static constexpr std::string_view extension = ".vtu";
    if (!vtuFilename.ends_with(extension))
    {
        vtuFilename += extension;
    }

    std::string vtpFilename = vtuFilename;
    const size_t dotPos = vtpFilename.rfind(".vtu");

    if (dotPos != std::string::npos)
    {
        vtpFilename.replace(dotPos, 4, "_boundary.vtp");
    }
    else
    {
        vtpFilename += "_boundary.vtp";
    }

    if (config.debug)
    {
        std::cout
            << '\n' << "Exporting results to VTK UnstructuredGrid..." << '\n';
    }

    VTK::writeVtkUnstructuredGrid
    (
        vtuFilename,
        mesh,
        scalarFieldsToVtk,
        vectorFieldsToVtk,
        config.debug
    );

    std::map<std::string, const FaceData<Scalar>*>
    boundaryScalarFields;

    for (const auto& output : turbulence.boundaryDataOutputs())
    {
        if (output.second != nullptr)
        {
            boundaryScalarFields[std::string{output.first}] =
                output.second;
        }
    }

    VTK::writeBoundaryData
    (
        vtpFilename,
        mesh,
        boundaryScalarFields,
        config.debug
    );

    std::cout
        << '\n' << "=== CFD Results Exported Successfully ===" << '\n';
    std::cout
        << "File: " << vtuFilename << '\n';
    std::cout
        << "Boundary: " << vtpFilename << '\n';
}

} // namespace PostProcess
