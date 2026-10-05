/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Forces.cpp
 * @brief Aerodynamic force calculation on a wall patch
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Forces.h"

// Standard library headers
#include <array>
#include <format>
#include <fstream>
#include <iostream>

// Project headers
#include "Vector.h"
#include "BoundaryConditions.h"
#include "BoundaryPatch.h"
#include "CaseConfiguration.h"
#include "ErrorHandler.h"
#include "Face.h"
#include "Field.h"
#include "Logger.h"
#include "Mesh.h"
#include "MomentumTransport.h"
#include "Comm.h"
#include "Reduce.h"
#include "TurbulenceModel.h"

// ***************************** namespace Forces ****************************

namespace Forces
{

// ***************************** Internal Helpers *****************************

namespace
{

// Drag/lift loads split into pressure and friction contributions
struct AeroForces
{
    Scalar pressureDrag;
    Scalar frictionDrag;
    Scalar pressureLift;
    Scalar frictionLift;
};

// Derive a forces output-file path from the configured output file name
FilePath forcesFilePath
(
    const FilePath& outputFilename,
    const FilePath& suffix
)
{
    static const std::array<FilePath, 4> extensions =
    {
        FilePath{".vtkhdf"},
        FilePath{".vtu"},
        FilePath{".vtp"},
        FilePath{".pvd"}
    };

    FilePath path = outputFilename;

    for (const FilePath& extension : extensions)
    {
        if (path.ends_with(extension))
        {
            path.resize(path.size() - extension.size());
            break;
        }
    }

    return path + suffix;
}

// Dynamic load 0.5 * rho * Vref^2 * A used to non-dimensionalize forces
Scalar referenceDynamicLoad(const CaseConfiguration& config)
{
    const Scalar referenceVelocity = magnitude(config.referenceVelocity);
    const Scalar dynamicLoad =
        S(0.5) * config.rho * referenceVelocity * referenceVelocity
      * config.referenceArea;

    if (dynamicLoad <= vSmallValue)
    {
        FatalError
        (
            "Reference dynamic load is zero; force coefficients are "
            "undefined. Check that density, referenceVelocity, and "
            "referenceArea are all non-zero."
        );
    }

    return dynamicLoad;
}

AeroForces computeForces
(
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    // Absent on this rank is fine; absent on every rank is a case error
    const BoundaryPatch* patch = nullptr;

    for (const BoundaryPatch& meshPatch : mesh.patches())
    {
        if (meshPatch.name() == config.forcesPatch)
        {
            patch = &meshPatch;
            break;
        }
    }

    if (!globalOr(patch != nullptr))
    {
        FatalError
        (
            "Forces patch '" + config.forcesPatch
          + "' not found on any rank"
        );
    }

    const ScalarField& Ux = solver.Ux();
    const ScalarField& Uy = solver.Uy();
    const ScalarField& Uz = solver.Uz();
    const ScalarField& pressure = solver.pressure();
    const FaceData<Scalar> wallShearStress =
        turbulence.wallShearStress(Ux, Uy, Uz);

    const FaceList& faces = mesh.faces();

    Vector pressureForceSum(S(0.0), S(0.0), S(0.0));
    Vector frictionForceSum(S(0.0), S(0.0), S(0.0));

    // A rank without the patch only joins the collective sums below
    if (patch != nullptr)
    {
        const Index firstFaceIdx = patch->firstFaceIdx();
        const Index lastFaceIdx = patch->lastFaceIdx();

        for
        (
            Index faceIdx = firstFaceIdx;
            faceIdx <= lastFaceIdx;
            ++faceIdx
        )
        {
            const Face& f = faces[faceIdx];
            const Index cellIdx = f.ownerCell();
            const Vector& normal = f.normal();

            // Pressure force from the kinematic pressure
            const Scalar pressureFace =
                bcManager.faceValue(f, pressure[cellIdx], Field::p);
            const Vector pressureContribution =
                (config.rho * pressureFace * f.projectedArea()) * normal;
            pressureForceSum += pressureContribution;

            // Skin-friction force
            const Vector cellVelocity(Ux[cellIdx], Uy[cellIdx], Uz[cellIdx]);
            const Scalar normalVelocity = dot(cellVelocity, normal);
            const Vector tangentVelocity =
                cellVelocity - normalVelocity * normal;
            const Scalar tangentMagnitude = magnitude(tangentVelocity);

            if (tangentMagnitude > vSmallValue)
            {
                const Vector shearDirection =
                    tangentVelocity / tangentMagnitude;
                const Scalar shearStress =
                    config.rho * wallShearStress[f.idx()];

                const Vector frictionContribution =
                    (shearStress * f.contactArea()) * shearDirection;
                frictionForceSum += frictionContribution;
            }
        }
    }

    // Each rank integrates its share of the patch; combine onto every rank
    const Vector pressureForce = globalSum(pressureForceSum);
    const Vector frictionForce = globalSum(frictionForceSum);

    const Vector& dragDir = config.dragDirection;
    const Vector& liftDir = config.liftDirection;

    return AeroForces
    {
        .pressureDrag = dot(pressureForce, dragDir),
        .frictionDrag = dot(frictionForce, dragDir),
        .pressureLift = dot(pressureForce, liftDir),
        .frictionLift = dot(frictionForce, liftDir)
    };
}

} // namespace

// ***************************** Force Reporting ******************************

void reportForces
(
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    const AeroForces forces =
        computeForces(solver, turbulence, mesh, bcManager, config);

    const Vector& dragDir = config.dragDirection;
    const Vector& liftDir = config.liftDirection;

    const Scalar pressureDrag = forces.pressureDrag;
    const Scalar frictionDrag = forces.frictionDrag;
    const Scalar totalDrag = pressureDrag + frictionDrag;

    const Scalar pressureLift = forces.pressureLift;
    const Scalar frictionLift = forces.frictionLift;
    const Scalar totalLift = pressureLift + frictionLift;

    // Non-dimensionalise: C = F / (0.5 * rho * Vref^2 * A).
    const Scalar dynamicLoad = referenceDynamicLoad(config);

    const Scalar pressureCd = pressureDrag / dynamicLoad;
    const Scalar frictionCd = frictionDrag / dynamicLoad;
    const Scalar totalCd = totalDrag / dynamicLoad;

    const Scalar pressureCl = pressureLift / dynamicLoad;
    const Scalar frictionCl = frictionLift / dynamicLoad;
    const Scalar totalCl = totalLift / dynamicLoad;

    // Write the breakdown to a text file beside the VTK output
    const FilePath outputPath =
        forcesFilePath(config.vtkOutputFilename, "_forces.txt");

    // Every rank reduced the forces above; only the master writes
    if (Comm::master())
    {
        std::ofstream file(outputPath);
        if (!file.is_open())
        {
            FatalError("Failed to open forces output file: " + outputPath);
        }

        file << std::format
        (
            "Aerodynamic forces\n"
            "Patch          : {}\n"
            "Drag direction : ({:.6f}, {:.6f}, {:.6f})\n"
            "Lift direction : ({:.6f}, {:.6f}, {:.6f})\n"
            "Reference U    : ({:.6f}, {:.6f}, {:.6f})\n"
            "\n"
            "{:<12}{:>16}{:>16}{:>16}\n"
            "{:<12}{:>16.6e}{:>16.6e}{:>16.6e}\n"
            "{:<12}{:>16.6e}{:>16.6e}{:>16.6e}\n"
            "\n"
            "Force coefficients (dimensionless), referenceArea = {:.6e}\n"
            "{:<12}{:>16}{:>16}{:>16}\n"
            "{:<12}{:>16.6e}{:>16.6e}{:>16.6e}\n"
            "{:<12}{:>16.6e}{:>16.6e}{:>16.6e}\n",
            config.forcesPatch,
            dragDir.x(), dragDir.y(), dragDir.z(),
            liftDir.x(), liftDir.y(), liftDir.z(),
            config.referenceVelocity.x(),
            config.referenceVelocity.y(),
            config.referenceVelocity.z(),
            "Force", "Pressure", "Friction", "Total",
            "Drag [N]", pressureDrag, frictionDrag, totalDrag,
            "Lift [N]", pressureLift, frictionLift, totalLift,
            config.referenceArea,
            "Coeff", "Pressure", "Friction", "Total",
            "Cd", pressureCd, frictionCd, totalCd,
            "Cl", pressureCl, frictionCl, totalCl
        );
        file.close();
    }

    // Print the breakdown to the console as two compact tables
    std::cout << '\n';
    Logger::sectionHeader("Aerodynamic Forces");
    Logger::subsection(std::format("Patch: {}", config.forcesPatch));

    Logger::breakdownHeader("Forces [N]");
    Logger::breakdownRow("Drag", pressureDrag, frictionDrag, totalDrag);
    Logger::breakdownRow("Lift", pressureLift, frictionLift, totalLift);

    Logger::breakdownHeader("Coefficients");
    Logger::breakdownRow("Cd", pressureCd, frictionCd, totalCd);
    Logger::breakdownRow("Cl", pressureCl, frictionCl, totalCl);

    Logger::subsection(std::format("Output file: {}", outputPath));
    Logger::iterationFooter();
}

// ************************* Transient Force History *************************

void writeForceHistoryHeader(const CaseConfiguration& config)
{
    const FilePath csvPath =
        forcesFilePath(config.vtkOutputFilename, "_forces.csv");

    if (!Comm::master())
    {
        return;
    }

    std::ofstream file(csvPath);
    if (!file.is_open())
    {
        FatalError("Failed to open forces history file: " + csvPath);
    }

    file << "time,pressureDrag,frictionDrag,totalDrag,"
         << "pressureLift,frictionLift,totalLift,Cd,Cl" << '\n';
}


void appendForceHistory
(
    Scalar time,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const CaseConfiguration& config
)
{
    const AeroForces forces =
        computeForces(solver, turbulence, mesh, bcManager, config);

    const Scalar totalDrag = forces.pressureDrag + forces.frictionDrag;
    const Scalar totalLift = forces.pressureLift + forces.frictionLift;

    const Scalar dynamicLoad = referenceDynamicLoad(config);

    const Scalar totalCd = totalDrag / dynamicLoad;
    const Scalar totalCl = totalLift / dynamicLoad;

    // Every rank reduced the forces above; only the master writes
    if (!Comm::master())
    {
        return;
    }

    const FilePath csvPath =
        forcesFilePath(config.vtkOutputFilename, "_forces.csv");

    std::ofstream file(csvPath, std::ios::app);
    if (!file.is_open())
    {
        FatalError("Failed to open forces history file: " + csvPath);
    }

    file << std::format
    (
        "{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n",
        time,
        forces.pressureDrag,
        forces.frictionDrag,
        totalDrag,
        forces.pressureLift,
        forces.frictionLift,
        totalLift,
        totalCd,
        totalCl
    );
}

} // namespace Forces
