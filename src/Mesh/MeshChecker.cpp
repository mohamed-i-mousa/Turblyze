/******************************************************************************
 * @file MeshChecker.cpp
 * @brief Mesh quality assessment utilities
 *****************************************************************************/

#include "MeshChecker.hpp"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numbers>

#include "ErrorHandler.hpp"


// ****************************** Public Methods ******************************

void MeshChecker::check() const
{
    std::cout
        << "\n--- Mesh Quality Check ---" << std::endl;

    if (mesh_.faces().empty() || mesh_.cells().empty())
    {
        FatalError("Empty mesh detected");
    }

    if (!validateConnectivity())
    {
        FatalError("Mesh connectivity validation failed");
    }

    // Face area statistics
    const Scalar firstArea = mesh_.faces()[0].projectedArea();
    Scalar minFaceArea = firstArea;
    Scalar maxFaceArea = firstArea;
    size_t minFaceIdx = mesh_.faces()[0].idx();
    size_t maxFaceIdx = mesh_.faces()[0].idx();

    // Collect faces with small area
    std::vector<size_t> smallAreaFaces;

    // Non-orthogonality statistics
    Scalar maxNonOrthogonality = 0.0;
    Scalar totalCosAngle = 0.0;
    size_t nonOrthCount = 0;
    size_t maxNonOrthFaceIdx = 0;
    std::vector<size_t> severeNonOrthFaces;

    // Skewness statistics
    Scalar maxSkewness = 0.0;
    size_t maxSkewFaceIdx = 0;
    std::vector<size_t> highSkewFaces;

    // Radian to degree conversion
    constexpr Scalar radToDeg = S(180.0) / std::numbers::pi_v<Scalar>;

    for (const auto& face : mesh_.faces())
    {
        const Scalar area = face.projectedArea();
        const size_t faceIdx = face.idx();

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

        if (area < minArea_)
        {
            smallAreaFaces.push_back(faceIdx);
        }

        // Calculate non-orthogonality and skewness
        if (face.isBoundary())
        {
            // Boundary face: calculate skewness only
            const Cell& ownerCell = mesh_.cells()[face.ownerCell()];

            Scalar skew = 
                calculateBoundarySkewness(face, ownerCell.centroid());

            if (skew > maxSkewness)
            {
                maxSkewness = skew;
                maxSkewFaceIdx = faceIdx;
            }
            if (skew > maxSkewThreshold_)
            {
                highSkewFaces.push_back(faceIdx);
            }
        }
        else
        {
            // Internal face: calculate both
            const Cell& ownerCell = mesh_.cells()[face.ownerCell()];

            const Cell& neighborCell =
                mesh_.cells()[face.neighborCell().value()];

            // Non-orthogonality (angle in degrees)
            Scalar ortho =
                calculateFaceOrthogonality
                (
                    ownerCell.centroid(),
                    neighborCell.centroid(),
                    face.normal()
                );

            Scalar angleRad = std::acos(ortho);
            Scalar angleDeg = angleRad * radToDeg;

            totalCosAngle += ortho;
            nonOrthCount++;

            if (angleDeg > maxNonOrthogonality)
            {
                maxNonOrthogonality = angleDeg;
                maxNonOrthFaceIdx = faceIdx;
            }

            if (angleDeg > maxNonOrthThreshold_)
            {
                severeNonOrthFaces.push_back(faceIdx);
            }

            // Skewness
            Scalar skew =
                calculateFaceSkewness
                (
                    face,
                    ownerCell.centroid(),
                    neighborCell.centroid()
                );

            if (skew > maxSkewness)
            {
                maxSkewness = skew;
                maxSkewFaceIdx = faceIdx;
            }

            if (skew > maxSkewThreshold_)
            {
                highSkewFaces.push_back(faceIdx);
            }
        }
    }

    // Calculate average non-orthogonality
    Scalar avgNonOrthogonality = std::acos(totalCosAngle / S(nonOrthCount));

    avgNonOrthogonality *= radToDeg;

    // Cell volume and aspect ratio statistics
    Scalar minCellVolume = mesh_.cells()[0].volume();
    Scalar maxCellVolume = mesh_.cells()[0].volume();
    size_t minCellIdx = mesh_.cells()[0].idx();
    size_t maxCellIdx = mesh_.cells()[0].idx();

    Scalar maxAspectRatio = 0.0;
    size_t maxAspectCellIdx = 0;
    std::vector<size_t> highAspectCells;

    std::vector<size_t> smallVolumeCells;
    std::vector<size_t> invertedCells;

    for (const auto& cell : mesh_.cells())
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
        else if (cell.volume() < minVolume_)
        {
            smallVolumeCells.push_back(cell.idx());
        }

        // Calculate aspect ratio
        Scalar aspectRatio = calculateCellAspectRatio(cell);
        if (aspectRatio > maxAspectRatio)
        {
            maxAspectRatio = aspectRatio;
            maxAspectCellIdx = cell.idx();
        }

        // High aspect ratio threshold
        if (aspectRatio > maxAspectThreshold_)
        {
            highAspectCells.push_back(cell.idx());
        }
    }

    // Store current format flags and precision
    const auto oldFlags = std::cout.flags();
    const auto oldPrecision = std::cout.precision();

    std::cout
        << "\nFace Area Statistics:" << std::endl;

    std::cout
        << "  Minimum area: "
        << std::scientific << std::setprecision(6)
        << minFaceArea << " m² (face "
        << minFaceIdx << ')' << std::endl;

    std::cout
        << "  Maximum area: "
        << std::scientific << std::setprecision(6)
        << maxFaceArea << " m² (face "
        << maxFaceIdx << ')' << std::endl;

    std::cout
        << "\nCell Volume Statistics:" << std::endl;

    std::cout
        << "  Minimum volume: "
        << std::scientific << std::setprecision(6)
        << minCellVolume << " m³ (cell "
        << minCellIdx << ')' << std::endl;

    std::cout
        << "  Maximum volume: "
        << std::scientific << std::setprecision(6)
        << maxCellVolume << " m³ (cell "
        << maxCellIdx << ')' << std::endl;

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
        << "\nNon-Orthogonality Statistics:" << std::endl;

    std::cout
        << std::fixed << std::setprecision(2);

    std::cout
        << "  Maximum: " << maxNonOrthogonality
        << "° (face " << maxNonOrthFaceIdx << ")"
        << std::endl;

    std::cout
        << "  Average: " << avgNonOrthogonality
        << "°" << std::endl;

    if (!severeNonOrthFaces.empty())
    {
        Warning
        (
            std::to_string(severeNonOrthFaces.size())
          + " faces with non-orthogonality > "
          + std::to_string(maxNonOrthThreshold_) + "°"
        );

        printIndicesList(severeNonOrthFaces, "Face");
    }

    // Skewness statistics
    std::cout
        << "\nSkewness Statistics:" << std::endl;

    std::cout
        << std::fixed << std::setprecision(3);

    std::cout
        << "  Maximum: " << maxSkewness << " (face "
        << maxSkewFaceIdx << ')' << std::endl;

    if (!highSkewFaces.empty())
    {
        Warning
        (
            std::to_string(highSkewFaces.size())
          + " faces with skewness > "
          + std::to_string(maxSkewThreshold_)
        );

        printIndicesList(highSkewFaces, "Face");
    }

    // Aspect ratio statistics
    std::cout
        << "\nAspect Ratio Statistics:" << std::endl;

    std::cout
        << std::fixed << std::setprecision(1);

    std::cout
        << "  Maximum: " << maxAspectRatio << " (cell "
        << maxAspectCellIdx << ')' << std::endl;

    if (!highAspectCells.empty())
    {
        Warning
        (
            std::to_string(highAspectCells.size())
          + " cells with aspect ratio > "
          + std::to_string(maxAspectThreshold_)
        );

        printIndicesList(highAspectCells, "Cell");
    }

    // Quality warnings for small areas/volumes
    if (!smallAreaFaces.empty())
    {
        std::cout
            << "\nQuality Check - Small Face Areas:"
            << std::endl;

        std::cout
            << "  Found " << smallAreaFaces.size()
            << " faces with area < ";

        std::cout
            << std::scientific
            << std::setprecision(0);

        std::cout
            << minArea_ << " m²" << std::endl;

        printIndicesList(smallAreaFaces, "Face");
    }

    if (!smallVolumeCells.empty())
    {
        std::cout
            << "\nQuality Check - Small Cell Volumes:"
            << std::endl;

        std::cout
            << "  Found " << smallVolumeCells.size()
            << " cells with volume < ";

        std::cout
            << std::scientific << std::setprecision(0);

        std::cout
            << minVolume_ << " m³" << std::endl;

        printIndicesList(smallVolumeCells, "Cell");
    }

    // Overall mesh quality summary
    std::cout
        << "\n--- Mesh Quality Summary ---" << std::endl;

    bool goodQuality = true;

    if (!invertedCells.empty())
    {
        Warning("Inverted cells detected");
        goodQuality = false;
    }

    if (maxNonOrthogonality > maxNonOrthThreshold_)
    {
        Warning
        (
            "Non-orthogonality exceeds "
          + std::to_string(maxNonOrthThreshold_)
          + "° threshold"
        );
        goodQuality = false;
    }

    if (maxSkewness > maxSkewThreshold_)
    {
        Warning
        (
            "Skewness exceeds "
          + std::to_string(maxSkewThreshold_)
          + " threshold"
        );
        goodQuality = false;
    }

    if (maxAspectRatio > maxAspectThreshold_)
    {
        Warning
        (
            "Aspect ratio exceeds "
          + std::to_string(maxAspectThreshold_)
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
            << " acceptable ranges"
            << std::endl;
    }

    // Restore original format flags and precision
    std::cout.flags(oldFlags);
    std::cout.precision(oldPrecision);
}


// ***************************** Private Methods ******************************

// Print up to 10 IDs from a list, with truncation indicator
void MeshChecker::printIndicesList
(
    std::span<const size_t> indices,
    std::string_view entityName
)
{
    constexpr size_t maxDisplay = 10;
    const size_t count =
        std::min(indices.size(), maxDisplay);

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

    for (size_t i = 0; i < count; ++i)
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
        << std::endl;
}


Scalar MeshChecker::calculateFaceOrthogonality
(
    const Vector& ownerCellCentroid,
    const Vector& neighborCellCentroid,
    const Vector& faceNormal
) const noexcept
{
    const Vector dPN = neighborCellCentroid - ownerCellCentroid;
    const Scalar dPNMag = dPN.magnitude();

    const Scalar cosAngle =
        dot(dPN, faceNormal) / (dPNMag + vSmallValue);

    return std::clamp(cosAngle, S(-1.0), S(1.0));
}


Scalar MeshChecker::calculateFaceSkewness
(
    const Face& face,
    const Vector& ownerCellCentroid,
    const Vector& neighborCellCentroid
) const
{
    const Vector dPf = face.centroid() - ownerCellCentroid;
    const Vector dPN = neighborCellCentroid - ownerCellCentroid;

    const Vector skewnessVector =
        dPf
      - ((dot(face.normal(), dPf))
      / (dot(face.normal(), dPN) + vSmallValue))
      * dPN;

    const Scalar skewnessMag = skewnessVector.magnitude();

    const Vector skewnessDirection =
        skewnessVector / (skewnessMag + smallValue);

    // Characteristic face dimension: empirical approximation
    Scalar faceCharacteristicLength = S(0.2) * dPN.magnitude() + vSmallValue;

    // Refine by finding max vertex extent in skewness direction
    std::span<const size_t> nodeIndices = face.nodeIndices();

    for (size_t nodeIdx : nodeIndices)
    {
        Vector vertexToCentroid = mesh_.nodes()[nodeIdx] - face.centroid();

        Scalar projection =
            std::abs(dot(skewnessDirection, vertexToCentroid));

        faceCharacteristicLength =
            std::max(faceCharacteristicLength, projection);
    }

    // Return normalized skewness
    return skewnessMag / faceCharacteristicLength;
}


Scalar MeshChecker::calculateBoundarySkewness
(
    const Face& face,
    const Vector& ownerCellCentroid
) const
{
    const Vector dPf = face.centroid() - ownerCellCentroid;

    // Virtual dPN for boundary faces
    const Vector dPN = S(2.0) * dot(face.normal(), dPf) * face.normal();

    const Vector skewnessVector = dPf - dot(face.normal(), dPf) * face.normal();

    const Scalar skewnessMag = skewnessVector.magnitude();

    const Vector skewnessDirection =
        skewnessVector / (skewnessMag + smallValue);

    // Characteristic face dimension: empirical approximation
    Scalar faceCharacteristicLength = S(0.4) * dPN.magnitude() + vSmallValue;

    // Refine by finding max vertex extent in skewness direction
    std::span<const size_t> nodeIndices = face.nodeIndices();

    for (size_t nodeIdx : nodeIndices)
    {
        Vector vertexToCentroid = mesh_.nodes()[nodeIdx] - face.centroid();

        Scalar projection = std::abs(dot(skewnessDirection, vertexToCentroid));

        faceCharacteristicLength =
            std::max(faceCharacteristicLength, projection);
    }

    return skewnessMag / faceCharacteristicLength;
}


Scalar MeshChecker::calculateCellAspectRatio
(
    const Cell& cell
) const
{
    // Accumulate absolute face area components per direction
    Vector sumMagAreaComponents(0.0, 0.0, 0.0);

    const auto& faceIndices = cell.faceIndices();

    for (size_t faceIdx : faceIndices)
    {
        const Face& face = mesh_.faces()[faceIdx];
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

bool MeshChecker::validateConnectivity() const
{
    bool valid = true;

    for (const auto& face : mesh_.faces())
    {
        if (face.ownerCell() >= mesh_.numCells())
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
         && face.neighborCell().value() >= mesh_.numCells())
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

        for (size_t nodeIdx : face.nodeIndices())
        {
            if (nodeIdx >= mesh_.numNodes())
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

    for (const auto& cell : mesh_.cells())
    {
        for (size_t faceIdx : cell.faceIndices())
        {
            if (faceIdx >= mesh_.numFaces())
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
