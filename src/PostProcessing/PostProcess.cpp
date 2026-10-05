/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
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
#include <limits>
#include <string>

// Project headers
#include "CaseConfiguration.h"
#include "Comm.h"
#include "DerivedFields.h"
#include "ErrorHandler.h"
#include "Logger.h"
#include "Mesh.h"
#include "MomentumTransport.h"
#include "Reduce.h"
#include "TurbulenceModel.h"
#include "HDF5BoundaryData.h"
#include "HDF5CellData.h"

// *************************** namespace PostProcess **************************

namespace PostProcess
{

// ***************************** Internal Helpers *****************************

namespace
{

// Output base path with any known trailing extension removed
FilePath outputBase(const FilePath& fileName)
{
    static const std::array<FilePath, 4> extensions =
    {
        FilePath{".vtkhdf"},
        FilePath{".vtu"},
        FilePath{".vtp"},
        FilePath{".pvd"}
    };

    FilePath base = fileName;

    for (const FilePath& extension : extensions)
    {
        if (base.ends_with(extension))
        {
            base.resize(base.size() - extension.size());
            break;
        }
    }

    return base;
}

} // namespace

void reportStatistics(const MomentumTransport& solver)
{
    std::cout << '\n';
    Logger::sectionHeader("Post-Processing Results");

    const ScalarField& Ux = solver.Ux();
    const ScalarField& Uy = solver.Uy();
    const ScalarField& Uz = solver.Uz();
    const ScalarField& pressure = solver.pressure();

    // Emptiness is decided globally: a per-rank return would deadlock
    const Count numDomainCells = solver.mesh().numDomainCells();
    const Count totalCells = globalSum(numDomainCells);

    if (totalCells == 0)
    {
        if (Comm::master())
        {
            Warning("Solution fields are empty. Skipping statistics.");
        }

        return;
    }

    const ScalarField velocityMag =
        VTK::velocityMagnitude(solver.mesh(), Ux, Uy, Uz);

    Scalar maximumVelocity = S(0.0);
    Scalar averageVelocity = S(0.0);
    Scalar maximumPressure = std::numeric_limits<Scalar>::lowest();
    Scalar minimumPressure = std::numeric_limits<Scalar>::max();

    for (Index cellIdx = 0; cellIdx < numDomainCells; ++cellIdx)
    {
        const Scalar vmag = velocityMag[cellIdx];
        maximumVelocity = std::max(maximumVelocity, vmag);
        averageVelocity += vmag;

        maximumPressure = std::max(maximumPressure, pressure[cellIdx]);
        minimumPressure = std::min(minimumPressure, pressure[cellIdx]);
    }
    maximumVelocity = globalMax(maximumVelocity);
    averageVelocity = globalSum(averageVelocity) / S(totalCells);
    maximumPressure = globalMax(maximumPressure);
    minimumPressure = globalMin(minimumPressure);

    Logger::subsection("Flow statistics");
    Logger::keyValue("Max velocity", maximumVelocity, "m/s");
    Logger::keyValue("Average velocity", averageVelocity, "m/s");
    Logger::keyValue("Pressure min", minimumPressure, "Pa");
    Logger::keyValue("Pressure max", maximumPressure, "Pa");
    Logger::iterationFooter();
}


void exportResults
(
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const CaseConfiguration& config
)
{
    std::cout << '\n';
    Logger::sectionHeader("Exporting Results");

    const FilePath cellDataFile = cellDataPath(config);
    const FilePath boundaryFile = boundaryDataPath(config);

    VTK::HDF5CellData volumeWriter
    (
        cellDataFile,
        mesh,
        config.debug
    );
    VTK::HDF5BoundaryData boundaryWriter
    (
        boundaryFile,
        mesh,
        &solver.bcManager(),
        config.debug
    );

    volumeWriter.writeGeometry();
    boundaryWriter.writeGeometry();

    // A steady result is a single-step VTKHDF series at t = 0
    appendTimeStep(volumeWriter, boundaryWriter, S(0.0), solver, turbulence);

    volumeWriter.finalize();
    boundaryWriter.finalize();

    Logger::keyValue("Cell data", cellDataFile);
    Logger::keyValue("Boundary data", boundaryFile);
    Logger::iterationFooter();
}


FilePath cellDataPath(const CaseConfiguration& config)
{
    return outputBase(config.vtkOutputFilename) + ".vtkhdf";
}


FilePath boundaryDataPath(const CaseConfiguration& config)
{
    return outputBase(config.vtkOutputFilename) + "_boundary.vtkhdf";
}


void appendTimeStep
(
    VTK::HDF5CellData& volumeWriter,
    VTK::HDF5BoundaryData& boundaryWriter,
    Scalar time,
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence
)
{
    const ScalarField& Ux = solver.Ux();
    const ScalarField& Uy = solver.Uy();
    const ScalarField& Uz = solver.Uz();
    const ScalarField& pressure = solver.pressure();

    const ScalarField velocityMag =
        VTK::velocityMagnitude(solver.mesh(), Ux, Uy, Uz);

    VTK::ScalarFieldMap scalarFieldsToVtk;
    scalarFieldsToVtk["pressure"] = &pressure;
    scalarFieldsToVtk["velocityMagnitude"] = &velocityMag;

    for (const auto& output : turbulence.cellDataOutputs())
    {
        scalarFieldsToVtk[Name{output.first}] = output.second;
    }

    VTK::VectorFieldMap vectorFieldsToVtk;
    vectorFieldsToVtk["velocity"] = {&Ux, &Uy, &Uz};

    volumeWriter.appendTimeStep(time, scalarFieldsToVtk, vectorFieldsToVtk);

    VTK::FaceDataMap boundaryScalarFields;

    for (const auto& output : turbulence.boundaryDataOutputs())
    {
        boundaryScalarFields[Name{output.first}] = output.second;
    }

    const FaceData<Scalar> wallShearStress =
        turbulence.wallShearStress(Ux, Uy, Uz);
    boundaryScalarFields["wallShearStress"] = &wallShearStress;

    boundaryWriter.appendTimeStep(time, boundaryScalarFields);
}

} // namespace PostProcess
