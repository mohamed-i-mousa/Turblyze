/******************************************************************************
 * @file MeshChecker.cpp
 * @brief Mesh quality assessment utilities
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header 
#include "MeshChecker.h"

// Standard library headers
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

// Project headers
#include "ErrorHandler.h"
#include "Integer.h"
#include "StringTypes.h"

// *************************** namespace MeshChecker **************************

namespace MeshChecker
{

// **************************** Quality Thresholds ****************************

constexpr Scalar minArea = S(1e-12);
constexpr Scalar minVolume = S(1e-30);
constexpr Scalar maxNonOrthThreshold = S(70.0);
constexpr Scalar maxSkewThreshold = S(4.0);
constexpr Scalar maxAspectThreshold = S(100.0);

// ***************************** Internal Helpers *****************************

namespace
{

// Print up to 10 IDs from a list, with truncation indicator
void printIndicesList
(
    IndexListRef indices,
    MessageRef entityName
)
{
    constexpr Count maxDisplay = 10;
    const Count count = std::min(indices.size(), maxDisplay);

    if (indices.size() <= maxDisplay)
    {
        std::cout
            << "  " << entityName << " IDs: ";
    }
    else
    {
        std::cout
            << "  First " << maxDisplay << " "
            << entityName << " IDs: ";
    }

    for (Index i = 0; i < count; ++i)
    {
        if (i > 0)
        {
            std::cout
                << ", ";
        }

        std::cout
            << indices[i];
    }

    if (indices.size() > maxDisplay)
    {
        std::cout
            << " ...";
    }

    std::cout
        << '\n';
}


Scalar faceOrthogonality
(
    const Vector& ownerCellCentroid,
    const Vector& neighborCellCentroid,
    const Vector& faceNormal
) noexcept
{
    const Vector dPN = neighborCellCentroid - ownerCellCentroid;
    const Scalar dPNMag = magnitude(dPN);
    const Scalar cosAngle = dot(dPN, faceNormal) / (dPNMag + vSmallValue);

    return std::clamp(cosAngle, S(-1.0), S(1.0));
}


Scalar faceSkewness
(
    const Mesh& mesh,
    const Face& face,
    const Vector& ownerCellCentroid,
    const Vector& neighborCellCentroid
)
{
    const Vector dPf = face.centroid() - ownerCellCentroid;
    const Vector dPN = neighborCellCentroid - ownerCellCentroid;

    const Vector skewnessVector =
        dPf
      - ((dot(face.normal(), dPf))
      / (dot(face.normal(), dPN) + vSmallValue))
      * dPN;

    const Scalar skewnessMag = magnitude(skewnessVector);

    const Vector skewnessDirection =
        skewnessVector / (skewnessMag + smallValue);

    // Characteristic face dimension: empirical approximation
    Scalar faceCharacteristicLength = S(0.2) * magnitude(dPN) + vSmallValue;

    // Refine by finding max vertex extent in skewness direction
    const IndexListRef nodeIndices = face.nodeIndices();

    for (Index nodeIdx : nodeIndices)
    {
        const Vector vertexToCentroid =
            mesh.nodes()[nodeIdx] - face.centroid();

        const Scalar projection =
            std::abs(dot(skewnessDirection, vertexToCentroid));

        faceCharacteristicLength =
            std::max(faceCharacteristicLength, projection);
    }

    // Return normalized skewness
    return skewnessMag / faceCharacteristicLength;
}


Scalar boundaryFaceSkewness
(
    const Mesh& mesh,
    const Face& face,
    const Vector& ownerCellCentroid
)
{
    const Vector dPf = face.centroid() - ownerCellCentroid;

    // Virtual dPN for boundary faces
    const Vector dPN = S(2.0) * dot(face.normal(), dPf) * face.normal();

    const Vector skewnessVector =
        dPf - dot(face.normal(), dPf) * face.normal();

    const Scalar skewnessMag = magnitude(skewnessVector);

    const Vector skewnessDirection =
        skewnessVector / (skewnessMag + smallValue);

    // Characteristic face dimension: empirical approximation
    Scalar faceCharacteristicLength = S(0.4) * magnitude(dPN) + vSmallValue;

    // Refine by finding max vertex extent in skewness direction
    const IndexListRef nodeIndices = face.nodeIndices();

    for (Index nodeIdx : nodeIndices)
    {
        const Vector vertexToCentroid =
            mesh.nodes()[nodeIdx] - face.centroid();

        const Scalar projection =
            std::abs(dot(skewnessDirection, vertexToCentroid));

        faceCharacteristicLength =
            std::max(faceCharacteristicLength, projection);
    }

    return skewnessMag / faceCharacteristicLength;
}


Scalar cellAspectRatio
(
    const Mesh& mesh,
    const Cell& cell
)
{
    // Accumulate absolute face area components per direction
    Vector sumMagAreaComponents;

    const auto faceIndices = cell.faceIndices();

    for (Index faceIdx : faceIndices)
    {
        const Face& face = mesh.faces()[faceIdx];
        const Vector areaVec = face.normal() * face.projectedArea();

        // Absolute components (sign irrelevant)
        sumMagAreaComponents.setX
        (
            sumMagAreaComponents.x() + std::abs(areaVec.x())
        );

        sumMagAreaComponents.setY
        (
            sumMagAreaComponents.y() + std::abs(areaVec.y())
        );

        sumMagAreaComponents.setZ
        (
            sumMagAreaComponents.z() + std::abs(areaVec.z())
        );
    }

    // Find min and max projected areas
    const Scalar minComponent =
        std::min
        (
            {
                sumMagAreaComponents.x(),
                sumMagAreaComponents.y(),
                sumMagAreaComponents.z()
            }
        );

    const Scalar maxComponent =
        std::max
        (
            {
                sumMagAreaComponents.x(),
                sumMagAreaComponents.y(),
                sumMagAreaComponents.z()
            }
        );

    Scalar directionalAspect = maxComponent / (minComponent + vSmallValue);

    // Add hydraulic aspect ratio for 3D cells
    const Scalar totalSurfaceArea =
        sumMagAreaComponents.x()
      + sumMagAreaComponents.y()
      + sumMagAreaComponents.z();

    const Scalar volume = cell.volume();

    if (volume > vSmallValue)
    {
        // Hydraulic aspect ratio: (1/6) * A / V^(2/3)
        const Scalar hydraulicAspect =
            (S(1.0)/S(6.0)) * totalSurfaceArea
          / std::pow(volume, S(2.0)/S(3.0));

        directionalAspect = std::max(directionalAspect, hydraulicAspect);
    }

    return directionalAspect;
}


bool validateConnectivity(const Mesh& mesh)
{
    bool valid = true;

    for (const auto& face : mesh.faces())
    {
        if (face.ownerCell() >= mesh.numCells())
        {
            Warning
            (
                "Face " + std::to_string(face.idx())
              + " owner cell index "
              + std::to_string(face.ownerCell())
              + " out of range"
            );
            valid = false;
        }

        if
        (
            !face.isBoundary()
         && face.neighborCell().value() >= mesh.numCells())
        {
            Warning
            (
                "Face " + std::to_string(face.idx())
              + " neighbor cell index "
              + std::to_string(face.neighborCell().value())
              + " out of range"
            );
            valid = false;
        }

        for (Index nodeIdx : face.nodeIndices())
        {
            if (nodeIdx >= mesh.numNodes())
            {
                Warning
                (
                    "Face " + std::to_string(face.idx())
                  + " node index " + std::to_string(nodeIdx)
                  + " out of range"
                );
                valid = false;
            }
        }
    }

    for (const auto& cell : mesh.cells())
    {
        for (Index faceIdx : cell.faceIndices())
        {
            if (faceIdx >= mesh.numFaces())
            {
                Warning
                (
                    "Cell " + std::to_string(cell.idx())
                  + " face index " + std::to_string(faceIdx)
                  + " out of range"
                );
                valid = false;
            }
        }
    }

    return valid;
}

}

// ******************************** Mesh Check ********************************

void check(const Mesh& mesh)
{
    std::cout
        << "\n--- Mesh Quality Check ---" << '\n';

    if (mesh.faces().empty() || mesh.cells().empty())
    {
        FatalError("Empty mesh detected");
    }

    if (!validateConnectivity(mesh))
    {
        FatalError("Mesh connectivity validation failed");
    }

    // Face area statistics
    const Scalar firstArea = mesh.faces()[0].projectedArea();
    Scalar minFaceArea = firstArea;
    Scalar maxFaceArea = firstArea;
    Index minFaceIdx = mesh.faces()[0].idx();
    Index maxFaceIdx = mesh.faces()[0].idx();

    // Collect faces with small area
    IndexList smallAreaFaces;

    // Non-orthogonality statistics
    Scalar maxNonOrthogonality = S(0.0);
    Scalar totalCosAngle = S(0.0);
    Count nonOrthCount = 0;
    Index maxNonOrthFaceIdx = 0;
    IndexList severeNonOrthFaces;

    // Skewness statistics
    Scalar maxSkewness = S(0.0);
    Index maxSkewFaceIdx = 0;
    IndexList highSkewFaces;

    // Radian to degree conversion
    constexpr Scalar radToDeg = S(180.0) / std::numbers::pi_v<Scalar>;

    for (const auto& face : mesh.faces())
    {
        const Scalar area = face.projectedArea();
        const Index faceIdx = face.idx();

        // Area statistics
        if (area < minFaceArea)
        {
            minFaceArea = area;
            minFaceIdx = faceIdx;
        }

        if (area > maxFaceArea)
        {
            maxFaceArea = area;
            maxFaceIdx = faceIdx;
        }

        if (area < minArea)
        {
            smallAreaFaces.push_back(faceIdx);
        }

        // Calculate non-orthogonality and skewness
        if (face.isBoundary())
        {
            // Boundary face: calculate skewness only
            const Cell& ownerCell = mesh.cells()[face.ownerCell()];

            const Scalar skew =
                boundaryFaceSkewness(mesh, face, ownerCell.centroid());

            if (skew > maxSkewness)
            {
                maxSkewness = skew;
                maxSkewFaceIdx = faceIdx;
            }
            if (skew > maxSkewThreshold)
            {
                highSkewFaces.push_back(faceIdx);
            }
        }
        else
        {
            // Internal face: calculate both
            const Cell& ownerCell = mesh.cells()[face.ownerCell()];

            const Cell& neighborCell =
                mesh.cells()[face.neighborCell().value()];

            // Non-orthogonality (angle in degrees)
            const Scalar ortho =
                faceOrthogonality
                (
                    ownerCell.centroid(),
                    neighborCell.centroid(),
                    face.normal()
                );

            const Scalar angleRad = std::acos(ortho);
            const Scalar angleDeg = angleRad * radToDeg;

            totalCosAngle += ortho;
            nonOrthCount++;

            if (angleDeg > maxNonOrthogonality)
            {
                maxNonOrthogonality = angleDeg;
                maxNonOrthFaceIdx = faceIdx;
            }

            if (angleDeg > maxNonOrthThreshold)
            {
                severeNonOrthFaces.push_back(faceIdx);
            }

            // Skewness
            const Scalar skew =
                faceSkewness
                (
                    mesh,
                    face,
                    ownerCell.centroid(),
                    neighborCell.centroid()
                );

            if (skew > maxSkewness)
            {
                maxSkewness = skew;
                maxSkewFaceIdx = faceIdx;
            }

            if (skew > maxSkewThreshold)
            {
                highSkewFaces.push_back(faceIdx);
            }
        }
    }

    // Calculate average non-orthogonality
    Scalar avgNonOrthogonality = S(0.0);
    if (nonOrthCount > 0)
    {
        avgNonOrthogonality =
            std::acos(totalCosAngle / S(nonOrthCount)) * radToDeg;
    }

    // Cell volume and aspect ratio statistics
    Scalar minCellVolume = mesh.cells()[0].volume();
    Scalar maxCellVolume = mesh.cells()[0].volume();
    Index minCellIdx = mesh.cells()[0].idx();
    Index maxCellIdx = mesh.cells()[0].idx();

    Scalar maxAspectRatio = S(0.0);
    Index maxAspectCellIdx = 0;
    IndexList highAspectCells;

    IndexList smallVolumeCells;
    IndexList invertedCells;

    for (const auto& cell : mesh.cells())
    {
        // Volume statistics
        if (cell.volume() < minCellVolume)
        {
            minCellVolume = cell.volume();
            minCellIdx = cell.idx();
        }

        if (cell.volume() > maxCellVolume)
        {
            maxCellVolume = cell.volume();
            maxCellIdx = cell.idx();
        }

        // Check for inverted or small cells
        if (cell.volume() < S(0.0))
        {
            invertedCells.push_back(cell.idx());
        }
        else if (cell.volume() < minVolume)
        {
            smallVolumeCells.push_back(cell.idx());
        }

        // Calculate aspect ratio
        const Scalar aspectRatio = cellAspectRatio(mesh, cell);
        if (aspectRatio > maxAspectRatio)
        {
            maxAspectRatio = aspectRatio;
            maxAspectCellIdx = cell.idx();
        }

        // High aspect ratio threshold
        if (aspectRatio > maxAspectThreshold)
        {
            highAspectCells.push_back(cell.idx());
        }
    }

    // Store current format flags and precision
    const auto oldFlags = std::cout.flags();
    const auto oldPrecision = std::cout.precision();

    std::cout
        << '\n' << "Face Area Statistics:" << '\n';

    std::cout
        << "  Minimum area: "
        << std::scientific << std::setprecision(6)
        << minFaceArea << " m² (face "
        << minFaceIdx << ')' << '\n';

    std::cout
        << "  Maximum area: "
        << std::scientific << std::setprecision(6)
        << maxFaceArea << " m² (face "
        << maxFaceIdx << ')' << '\n';

    std::cout
        << '\n' << "Cell Volume Statistics:" << '\n';

    std::cout
        << "  Minimum volume: "
        << std::scientific << std::setprecision(6)
        << minCellVolume << " m³ (cell "
        << minCellIdx << ')' << '\n';

    std::cout
        << "  Maximum volume: "
        << std::scientific << std::setprecision(6)
        << maxCellVolume << " m³ (cell "
        << maxCellIdx << ')' << '\n';

    if (!invertedCells.empty())
    {
        FatalError
        (
            std::to_string(invertedCells.size())
          + " inverted cells (negative volume) detected"
        );
    }

    // Non-orthogonality statistics
    std::cout
        << '\n' << "Non-Orthogonality Statistics:" << '\n';

    std::cout
        << std::fixed << std::setprecision(2);

    if (nonOrthCount > 0)
    {
        std::cout
            << "  Maximum: " << maxNonOrthogonality
            << "° (face " << maxNonOrthFaceIdx << ")" << '\n';

        std::cout
            << "  Average: " << avgNonOrthogonality << "°" << '\n';
    }
    else
    {
        std::cout
            << "  No internal faces to measure" << '\n';
    }

    if (!severeNonOrthFaces.empty())
    {
        Warning
        (
            std::to_string(severeNonOrthFaces.size())
          + " faces with non-orthogonality > "
          + std::to_string(maxNonOrthThreshold) + "°"
        );

        printIndicesList(severeNonOrthFaces, "Face");
    }

    // Skewness statistics
    std::cout
        << '\n' << "Skewness Statistics:" << '\n';

    std::cout
        << std::fixed << std::setprecision(3);

    std::cout
        << "  Maximum: " << maxSkewness << " (face "
        << maxSkewFaceIdx << ')' << '\n';

    if (!highSkewFaces.empty())
    {
        Warning
        (
            std::to_string(highSkewFaces.size())
          + " faces with skewness > "
          + std::to_string(maxSkewThreshold)
        );

        printIndicesList(highSkewFaces, "Face");
    }

    // Aspect ratio statistics
    std::cout
        << '\n' << "Aspect Ratio Statistics:" << '\n';

    std::cout
        << std::fixed << std::setprecision(1);

    std::cout
        << "  Maximum: " << maxAspectRatio << " (cell "
        << maxAspectCellIdx << ')' << '\n';

    if (!highAspectCells.empty())
    {
        Warning
        (
            std::to_string(highAspectCells.size())
          + " cells with aspect ratio > "
          + std::to_string(maxAspectThreshold)
        );

        printIndicesList(highAspectCells, "Cell");
    }

    // Quality warnings for small areas/volumes
    if (!smallAreaFaces.empty())
    {
        std::cout
            << '\n' << "Quality Check - Small Face Areas:" << '\n';

        std::cout
            << "  Found " << smallAreaFaces.size() << " faces with area < ";

        std::cout
            << std::scientific
            << std::setprecision(0);

        std::cout
            << minArea << " m²" << '\n';

        printIndicesList(smallAreaFaces, "Face");
    }

    if (!smallVolumeCells.empty())
    {
        std::cout
            << '\n' << "Quality Check - Small Cell Volumes:" << '\n';

        std::cout
            << "  Found " << smallVolumeCells.size()
            << " cells with volume < ";

        std::cout
            << std::scientific << std::setprecision(0);

        std::cout
            << minVolume << " m³" << '\n';

        printIndicesList(smallVolumeCells, "Cell");
    }

    // Overall mesh quality summary
    std::cout
        << "\n--- Mesh Quality Summary ---" << '\n';

    bool goodQuality = true;

    if (!invertedCells.empty())
    {
        Warning("Inverted cells detected");
        goodQuality = false;
    }

    if (maxNonOrthogonality > maxNonOrthThreshold)
    {
        Warning
        (
            "Non-orthogonality exceeds "
          + std::to_string(maxNonOrthThreshold)
          + "° threshold"
        );
        goodQuality = false;
    }

    if (maxSkewness > maxSkewThreshold)
    {
        Warning
        (
            "Skewness exceeds "
          + std::to_string(maxSkewThreshold)
          + " threshold"
        );
        goodQuality = false;
    }

    if (maxAspectRatio > maxAspectThreshold)
    {
        Warning
        (
            "Aspect ratio exceeds "
          + std::to_string(maxAspectThreshold)
          + " threshold"
        );
        goodQuality = false;
    }

    if
    (
        smallAreaFaces.empty()
     && smallVolumeCells.empty()
     && invertedCells.empty()
     && goodQuality
    )
    {
        std::cout
            << "DONE: All mesh quality metrics within"
            << " acceptable ranges" << '\n';
    }

    // Restore original format flags and precision
    std::cout.flags(oldFlags);
    std::cout.precision(oldPrecision);
}

} // namespace MeshChecker
