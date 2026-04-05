/******************************************************************************
 * @file MeshChecker.cpp
 * @brief Mesh quality assessment utilities
 *****************************************************************************/

#include "MeshChecker.hpp"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

#include "ErrorHandler.hpp"
#include <numbers>
#include <stdexcept>


// ************************* Constructor & Destructor *************************

MeshChecker::MeshChecker
(
    std::span<const Vector> nodes,
    std::span<const Face> faces,
    std::span<const Cell> cells
) noexcept
:
    allNodes_(nodes),
    allFaces_(faces),
    allCells_(cells)
{}


// ****************************** Public Methods ******************************

void MeshChecker::check() const
{
    std::cout
        << "\n--- Mesh Quality Check ---" << std::endl;

    if (allFaces_.empty() || allCells_.empty())
    {
        std::cout
            << "Warning: Empty mesh detected!" << std::endl;
        return;
    }

    if (!validateConnectivity())
    {
        std::cout
            << "Mesh connectivity validation failed."
            << " Skipping quality checks."
            << std::endl;
        return;
    }

    // Face area statistics
    const Scalar firstArea = allFaces_[0].projectedArea();
    Scalar minFaceArea = firstArea;
    Scalar maxFaceArea = firstArea;
    size_t minFaceId = allFaces_[0].idx();
    size_t maxFaceId = allFaces_[0].idx();

    std::vector<size_t> smallAreaFaces;

    // Non-orthogonality statistics
    Scalar maxNonOrthogonality = 0.0;
    Scalar totalCosAngle = 0.0;
    size_t nonOrthCount = 0;
    size_t maxNonOrthFaceId = 0;
    std::vector<size_t> severeNonOrthFaces;

    // Skewness statistics
    Scalar maxSkewness = 0.0;
    size_t maxSkewFaceId = 0;
    std::vector<size_t> highSkewFaces;

    // Radian to degree conversion
    constexpr Scalar radToDeg = S(180.0) / std::numbers::pi_v<Scalar>;

    for (const auto& face : allFaces_)
    {
        const Scalar area = face.projectedArea();
        const size_t faceIdx = face.idx();

        // Area statistics
        if (area < minFaceArea)
        {
            minFaceArea = area;
            minFaceId = faceIdx;
        }
        if (area > maxFaceArea)
        {
            maxFaceArea = area;
            maxFaceId = faceIdx;
        }

        if (area < minArea_)
        {
            smallAreaFaces.push_back(faceIdx);
        }

        // Calculate non-orthogonality and skewness
        if (face.isBoundary())
        {
            // Boundary face: calculate skewness only
            const Cell& ownerCell = allCells_[face.ownerCell()];

            Scalar skew =
                calculateBoundarySkewness
                (
                    face,
                    ownerCell.centroid(),
                    face.centroid(),
                    face.normal()
                );

            if (skew > maxSkewness)
            {
                maxSkewness = skew;
                maxSkewFaceId = faceIdx;
            }
            if (skew > maxSkewThreshold_)
            {
                highSkewFaces.push_back(faceIdx);
            }
        }
        else
        {
            // Internal face: calculate both
            const Cell& ownerCell =
                allCells_[face.ownerCell()];

            const Cell& neighborCell =
                allCells_[face.neighborCell().value()];

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
                maxNonOrthFaceId = faceIdx;
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
                    neighborCell.centroid(),
                    face.centroid(),
                    face.normal()
                );

            if (skew > maxSkewness)
            {
                maxSkewness = skew;
                maxSkewFaceId = faceIdx;
            }

            if (skew > maxSkewThreshold_)
            {
                highSkewFaces.push_back(faceIdx);
            }
        }
    }

    // Calculate average non-orthogonality
    Scalar avgNonOrthogonality =
        std::acos(totalCosAngle / S(nonOrthCount));

    avgNonOrthogonality *= radToDeg;

    // Cell volume and aspect ratio statistics
    Scalar minCellVolume = allCells_[0].volume();
    Scalar maxCellVolume = allCells_[0].volume();
    size_t minCellId = allCells_[0].idx();
    size_t maxCellId = allCells_[0].idx();

    Scalar maxAspectRatio = 0.0;
    size_t maxAspectCellId = 0;
    std::vector<size_t> highAspectCells;

    std::vector<size_t> smallVolumeCells;
    std::vector<size_t> invertedCells;

    for (const auto& cell : allCells_)
    {
        // Volume statistics
        if (cell.volume() < minCellVolume)
        {
            minCellVolume = cell.volume();
            minCellId = cell.idx();
        }
        if (cell.volume() > maxCellVolume)
        {
            maxCellVolume = cell.volume();
            maxCellId = cell.idx();
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
            maxAspectCellId = cell.idx();
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
        << minFaceId << ')' << std::endl;

    std::cout
        << "  Maximum area: "
        << std::scientific << std::setprecision(6)
        << maxFaceArea << " m² (face "
        << maxFaceId << ')' << std::endl;

    std::cout
        << "\nCell Volume Statistics:" << std::endl;

    std::cout
        << "  Minimum volume: "
        << std::scientific << std::setprecision(6)
        << minCellVolume << " m³ (cell "
        << minCellId << ')' << std::endl;

    std::cout
        << "  Maximum volume: "
        << std::scientific << std::setprecision(6)
        << maxCellVolume << " m³ (cell "
        << maxCellId << ')' << std::endl;

    if (!invertedCells.empty())
    {
        std::cout
            << "\nCRITICAL: "
            << invertedCells.size()
            << " inverted cells (negative volume)!"
            << std::endl;

        printIndicesList(invertedCells, "Cell");
    }

    // Non-orthogonality statistics
    std::cout
        << "\nNon-Orthogonality Statistics:" << std::endl;

    std::cout
        << std::fixed << std::setprecision(2);

    std::cout
        << "  Maximum: " << maxNonOrthogonality
        << "° (face " << maxNonOrthFaceId << ")"
        << std::endl;

    std::cout
        << "  Average: " << avgNonOrthogonality
        << "°" << std::endl;

    if (!severeNonOrthFaces.empty())
    {
        std::cout
            << "  WARNING: " << severeNonOrthFaces.size()
            << " faces with non-orthogonality > "
            << maxNonOrthThreshold_ << "°"
            << std::endl;

        printIndicesList(severeNonOrthFaces, "Face");
    }

    // Skewness statistics
    std::cout
        << "\nSkewness Statistics:" << std::endl;

    std::cout
        << std::fixed << std::setprecision(3);

    std::cout
        << "  Maximum: " << maxSkewness << " (face "
        << maxSkewFaceId << ')' << std::endl;

    if (!highSkewFaces.empty())
    {
        std::cout
            << "  WARNING: " << highSkewFaces.size()
            << " faces with skewness > "
            << maxSkewThreshold_
            << std::endl;

        printIndicesList(highSkewFaces, "Face");
    }

    // Aspect ratio statistics
    std::cout
        << "\nAspect Ratio Statistics:" << std::endl;

    std::cout
        << std::fixed << std::setprecision(1);

    std::cout
        << "  Maximum: " << maxAspectRatio << " (cell "
        << maxAspectCellId << ')' << std::endl;

    if (!highAspectCells.empty())
    {
        std::cout
            << "  WARNING: " << highAspectCells.size()
            << " cells with aspect ratio > "
            << maxAspectThreshold_
            << std::endl;

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
        std::cout
            << "WARNING: Inverted cells detected"
            << std::endl;

        goodQuality = false;
    }

    if (maxNonOrthogonality > maxNonOrthThreshold_)
    {
        std::cout
            << "WARNING: Non-orthogonality exceeds "
            << maxNonOrthThreshold_
            << "° threshold" << std::endl;

        goodQuality = false;
    }

    if (maxSkewness > maxSkewThreshold_)
    {
        std::cout
            << "WARNING: Skewness exceeds "
            << maxSkewThreshold_
            << " threshold" << std::endl;

        goodQuality = false;
    }

    if (maxAspectRatio > maxAspectThreshold_)
    {
        std::cout
            << "WARNING: Aspect ratio exceeds "
            << maxAspectThreshold_
            << " threshold" << std::endl;

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
    const Vector& neighborCellCentroid,
    const Vector& faceCentroid,
    const Vector& faceNormal
) const
{
    const Vector dPf = faceCentroid - ownerCellCentroid;
    const Vector dPN = neighborCellCentroid - ownerCellCentroid;

    const Vector skewnessVector =
        dPf
      - ((dot(faceNormal, dPf))
      / (dot(faceNormal, dPN) + vSmallValue))
      * dPN;

    const Scalar skewnessMag = skewnessVector.magnitude();

    const Vector skewnessDirection =
        skewnessVector / (skewnessMag + smallValue);

    // Characteristic face dimension: empirical approximation
    Scalar faceCharacteristicLength = S(0.2) * dPN.magnitude() + vSmallValue;

    // Refine by finding max vertex extent in skewness dir
    std::span<const size_t> nodeIndices = face.nodeIndices();

    for (size_t nodeIdx : nodeIndices)
    {
        if (nodeIdx >= allNodes_.size())
        {
            FatalError
            (
                "Node index "
              + std::to_string(nodeIdx)
              + " out of range in skewness"
                " calculation for face "
              + std::to_string(face.idx())
            );
        }

        Vector vertexToCentroid =
            allNodes_[nodeIdx] - faceCentroid;

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
    const Vector& ownerCellCentroid,
    const Vector& faceCentroid,
    const Vector& faceNormal
) const
{
    const Vector dPf = faceCentroid - ownerCellCentroid;

    // Virtual dPN for boundary faces
    const Vector dPN = S(2.0) * dot(faceNormal, dPf) * faceNormal;

    const Vector skewnessVector = dPf - dot(faceNormal, dPf) * faceNormal;

    const Scalar skewnessMag = skewnessVector.magnitude();

    const Vector skewnessDirection =
        skewnessVector / (skewnessMag + smallValue);

    // Characteristic face dimension: empirical approximation
    Scalar faceCharacteristicLength = S(0.4) * dPN.magnitude() + vSmallValue;

    // Refine by finding max vertex extent in skewness dir
    std::span<const size_t> nodeIndices = face.nodeIndices();

    for (size_t nodeIdx : nodeIndices)
    {
        if (nodeIdx >= allNodes_.size())
        {
            FatalError
            (
                "Node index "
              + std::to_string(nodeIdx)
              + " out of range in boundary"
                " skewness calculation for face "
              + std::to_string(face.idx())
            );
        }

        Vector vertexToCentroid = allNodes_[nodeIdx] - faceCentroid;

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
        const Face& face = allFaces_[faceIdx];
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

    for (const auto& face : allFaces_)
    {
        if (face.ownerCell() >= allCells_.size())
        {
            std::cout
                << "ERROR: Face " << face.idx()
                << " owner cell index "
                << face.ownerCell() << " out of range"
                << std::endl;
            valid = false;
        }

        if (!face.isBoundary()
         && face.neighborCell().value() >= allCells_.size())
        {
            std::cout
                << "ERROR: Face " << face.idx()
                << " neighbor cell index "
                << face.neighborCell().value()
                << " out of range"
                << std::endl;
            valid = false;
        }

        for (size_t nodeIdx : face.nodeIndices())
        {
            if (nodeIdx >= allNodes_.size())
            {
                std::cout
                    << "ERROR: Face " << face.idx()
                    << " node index " << nodeIdx
                    << " out of range"
                    << std::endl;
                valid = false;
            }
        }
    }

    for (const auto& cell : allCells_)
    {
        for (size_t faceIdx : cell.faceIndices())
        {
            if (faceIdx >= allFaces_.size())
            {
                std::cout
                    << "ERROR: Cell " << cell.idx()
                    << " face index " << faceIdx
                    << " out of range"
                    << std::endl;
                valid = false;
            }
        }
    }

    return valid;
}
