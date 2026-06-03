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

    using ScalarFieldMap = std::map<Name, const ScalarField*>;
    using VectorFieldMap = std::map<Name, std::array<const ScalarField*, 3>>;
    using BoundaryScalarFieldMap = std::map<Name, const FaceData<Scalar>*>;

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

    for (Index cellIdx = 0; cellIdx < Ux.size(); ++cellIdx)
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

    ScalarFieldMap scalarFieldsToVtk;

    scalarFieldsToVtk["pressure"] = &pressure;
    scalarFieldsToVtk["velocityMagnitude"] = &velocityMag;

    for (const auto& output : turbulence.cellDataOutputs())
    {
        if (output.second != nullptr)
        {
            scalarFieldsToVtk[Name{output.first}] = output.second;
        }
    }

    VectorFieldMap vectorFieldsToVtk;

    vectorFieldsToVtk["velocity"] = {&Ux, &Uy, &Uz};

    FilePath vtuFilename = config.vtkOutputFilename;
    static constexpr FilePathRef extension = ".vtu";
    if (!vtuFilename.ends_with(extension))
    {
        vtuFilename += extension;
    }

    FilePath vtpFilename = vtuFilename;
    const Index dotPos = vtpFilename.rfind(".vtu");

    if (dotPos != FilePath::npos)
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

    BoundaryScalarFieldMap boundaryScalarFields;

    for (const auto& output : turbulence.boundaryDataOutputs())
    {
        if (output.second != nullptr)
        {
            boundaryScalarFields[Name{output.first}] =
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
