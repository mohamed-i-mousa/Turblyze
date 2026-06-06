/******************************************************************************
 * @file Forces.cpp
 * @brief Aerodynamic force calculation on a wall patch
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Forces.h"

// Standard library headers
#include <fstream>
#include <iomanip>
#include <iostream>

// Project headers
#include "BoundaryConditions.h"
#include "BoundaryPatch.h"
#include "CaseConfiguration.h"
#include "ErrorHandler.h"
#include "Face.h"
#include "Field.h"
#include "Logger.h"
#include "Mesh.h"
#include "SIMPLE.h"
#include "TurbulenceModel.h"
#include "Vector.h"

// ***************************** namespace Forces ****************************

namespace Forces
{

// ***************************** Internal Helpers *****************************

namespace
{

// Derive the forces text-file path from the VTK output filename
FilePath forcesFilePath(const FilePath& vtkOutputFilename)
{
    FilePath path = vtkOutputFilename;
    const Index dotPos = path.rfind(".vtu");

    if (dotPos != FilePath::npos)
    {
        path.replace(dotPos, 4, "_forces.txt");
    }
    else
    {
        path += "_forces.txt";
    }

    return path;
}

} // namespace

// ***************************** Force Reporting ******************************

void reportForces
(
    const SIMPLE& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    const BoundaryPatch& patch = bcManager.patch(config.forcesPatch);

    const ScalarField& Ux = solver.Ux();
    const ScalarField& Uy = solver.Uy();
    const ScalarField& Uz = solver.Uz();
    const ScalarField& pressure = solver.pressure();
    const FaceData<Scalar> wallShearStress =
        turbulence.wallShearStress(Ux, Uy, Uz);

    const FaceListRef faces = mesh.faces();

    // Integrated dynamic force vectors [N] over the patch
    Vector pressureForce;
    Vector frictionForce;

    for
    (
        Index faceIdx = patch.firstFaceIdx();
        faceIdx <= patch.lastFaceIdx();
        ++faceIdx
    )
    {
        const Face& face = faces[faceIdx];
        const Index cellIdx = face.ownerCell();
        const Vector& normal = face.normal();

        // Pressure force from the kineamtic pressure
        const Scalar pressureFace =
            bcManager.boundaryFaceValue(Field::p, pressure, face);
        pressureForce +=
            (config.rho * pressureFace * face.projectedArea()) * normal;

        // Skin-friction force
        const Vector cellVelocity(Ux[cellIdx], Uy[cellIdx], Uz[cellIdx]);
        const Scalar normalVelocity = dot(cellVelocity, normal);
        const Vector tangentVelocity =
            cellVelocity - normalVelocity * normal;
        const Scalar tangentMagnitude = magnitude(tangentVelocity);

        if (tangentMagnitude > vSmallValue)
        {
            const Vector shearDirection = tangentVelocity / tangentMagnitude;
            const Scalar shearStress = config.rho * wallShearStress[face.idx()];

            frictionForce +=
                (shearStress * face.contactArea()) * shearDirection;
        }
    }

    // Decompose onto the drag and lift directions
    const Vector& dragDir = config.dragDirection;
    const Vector& liftDir = config.liftDirection;

    const Scalar pressureDrag = dot(pressureForce, dragDir);
    const Scalar frictionDrag = dot(frictionForce, dragDir);
    const Scalar totalDrag = pressureDrag + frictionDrag;

    const Scalar pressureLift = dot(pressureForce, liftDir);
    const Scalar frictionLift = dot(frictionForce, liftDir);
    const Scalar totalLift = pressureLift + frictionLift;

    // Non-dimensionalise: C = F / (0.5 * rho * Vref^2 * A).
    const Scalar referenceVelocity = magnitude(config.referenceVelocity);
    const Scalar dynamicLoad =
        S(0.5) * config.rho * referenceVelocity * referenceVelocity
      * config.referenceArea;

    const Scalar pressureCd = pressureDrag / dynamicLoad;
    const Scalar frictionCd = frictionDrag / dynamicLoad;
    const Scalar totalCd = totalDrag / dynamicLoad;

    const Scalar pressureCl = pressureLift / dynamicLoad;
    const Scalar frictionCl = frictionLift / dynamicLoad;
    const Scalar totalCl = totalLift / dynamicLoad;

    // Write the breakdown to a text file beside the VTK output
    const FilePath outputPath = forcesFilePath(config.vtkOutputFilename);

    std::ofstream file(outputPath);
    if (!file.is_open())
    {
        FatalError("Failed to open forces output file: " + outputPath);
    }

    file << std::scientific << std::setprecision(6);
    file
        << "Aerodynamic forces" << '\n'
        << "Patch          : " << config.forcesPatch << '\n'
        << "Drag direction : " << dragDir << '\n'
        << "Lift direction : " << liftDir << '\n'
        << "Reference U    : " << config.referenceVelocity << '\n'
        << '\n'
        << std::left << std::setw(12) << "Force"
        << std::right
        << std::setw(16) << "Pressure"
        << std::setw(16) << "Friction"
        << std::setw(16) << "Total" << '\n'
        << std::left << std::setw(12) << "Drag [N]"
        << std::right
        << std::setw(16) << pressureDrag
        << std::setw(16) << frictionDrag
        << std::setw(16) << totalDrag << '\n'
        << std::left << std::setw(12) << "Lift [N]"
        << std::right
        << std::setw(16) << pressureLift
        << std::setw(16) << frictionLift
        << std::setw(16) << totalLift << '\n'
        << '\n'
        << "Force coefficients (dimensionless), referenceArea = "
        << config.referenceArea << '\n'
        << std::left << std::setw(12) << "Coeff"
        << std::right
        << std::setw(16) << "Pressure"
        << std::setw(16) << "Friction"
        << std::setw(16) << "Total" << '\n'
        << std::left << std::setw(12) << "Cd"
        << std::right
        << std::setw(16) << pressureCd
        << std::setw(16) << frictionCd
        << std::setw(16) << totalCd << '\n'
        << std::left << std::setw(12) << "Cl"
        << std::right
        << std::setw(16) << pressureCl
        << std::setw(16) << frictionCl
        << std::setw(16) << totalCl << '\n';
    file.close();

    // Print the breakdown to the console as two compact tables
    std::cout << '\n';
    Logger::sectionHeader("Aerodynamic Forces");
    Logger::subsection("Patch: " + config.forcesPatch);

    Logger::breakdownHeader("Forces [N]");
    Logger::breakdownRow("Drag", pressureDrag, frictionDrag, totalDrag);
    Logger::breakdownRow("Lift", pressureLift, frictionLift, totalLift);

    Logger::breakdownHeader("Coefficients");
    Logger::breakdownRow("Cd", pressureCd, frictionCd, totalCd);
    Logger::breakdownRow("Cl", pressureCl, frictionCl, totalCl);

    Logger::subsection("Output file: " + outputPath);
    Logger::iterationFooter();
}

} // namespace Forces
