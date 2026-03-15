/******************************************************************************
 * @file MeshChecker.cpp
 * @brief Mesh quality assessment utilities
 *****************************************************************************/

#include "MeshChecker.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <stdexcept>


// ************************* Constructor & Destructor *************************

MeshChecker::MeshChecker
(
    const std::vector<Vector>& nodes,
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells
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
            if (face.ownerCell() >= allCells_.size())
            {
                throw
                    std::out_of_range
                    (
                        "Owner cell index "
                      + std::to_string(face.ownerCell())
                      + " out of range for face "
                      + std::to_string(faceIdx)
                    );
            }

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
            if (face.ownerCell() >= allCells_.size())
            {
                throw
                    std::out_of_range
                    (
                        "Owner cell index "
                      + std::to_string(face.ownerCell())
                      + " out of range for face "
                      + std::to_string(faceIdx)
                    );
            }

            if (face.neighborCell().value() >= allCells_.size())
            {
                throw
                    std::out_of_range
                    (
                        "Neighbor cell index "
                      + std::to_string(
                            face.neighborCell().value()
                        )
                      + " out of range for face "
                      + std::to_string(faceIdx)
                    );
            }

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
            Scalar angleDeg =
                angleRad * S(180.0) / S(M_PI);

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

    avgNonOrthogonality = avgNonOrthogonality * S(180.0) / S(M_PI);

    // Cell volume and aspect ratio statistics
    Scalar minCellVolume = allCells_[0].volume();
    Scalar maxCellVolume = allCells_[0].volume();
    size_t minCellId = allCells_[0].idx();
    size_t maxCellId = allCells_[0].idx();

    Scalar maxAspectRatio = 0.0;
    size_t maxAspectCellId = 0;
    std::vector<size_t> highAspectCells;

    std::vector<size_t> smallVolumeCells;

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

        // Check for cells with volume smaller than 1e-30
        if (cell.volume() < minVolume_)
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
        << minFaceId << ")" << std::endl;

    std::cout
        << "  Maximum area: "
        << std::scientific << std::setprecision(6)
        << maxFaceArea << " m² (face "
        << maxFaceId << ")" << std::endl;

    std::cout
        << "\nCell Volume Statistics:" << std::endl;

    std::cout
        << "  Minimum volume: "
        << std::scientific << std::setprecision(6)
        << minCellVolume << " m³ (cell "
        << minCellId << ")" << std::endl;

    std::cout
        << "  Maximum volume: "
        << std::scientific << std::setprecision(6)
        << maxCellVolume << " m³ (cell "
        << maxCellId << ")" << std::endl;

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
        << maxSkewFaceId << ")" << std::endl;

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
        << maxAspectCellId << ")" << std::endl;

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

    if (maxNonOrthogonality > maxNonOrthThreshold_)
    {
        std::cout
            << "⚠ Non-orthogonality exceeds "
            << maxNonOrthThreshold_
            << "° threshold" << std::endl;

        goodQuality = false;
    }

    if (maxSkewness > maxSkewThreshold_)
    {
        std::cout
            << "⚠ Skewness exceeds "
            << maxSkewThreshold_
            << " threshold" << std::endl;

        goodQuality = false;
    }

    if (maxAspectRatio > maxAspectThreshold_)
    {
        std::cout
            << "⚠ Aspect ratio exceeds "
            << maxAspectThreshold_
            << " threshold" << std::endl;

        goodQuality = false;
    }

    if (smallAreaFaces.empty()
     && smallVolumeCells.empty()
     && goodQuality)
    {
        std::cout
            << "✓ All mesh quality metrics within"
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
    const std::vector<size_t>& indices,
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
    const Scalar faceNormalMag = faceNormal.magnitude();

    const Scalar cosAngle =
        dot(dPN, faceNormal)
      / (dPNMag * faceNormalMag + vSmallValue);

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
    const std::vector<size_t>& nodeIndices = face.nodeIndices();

    for (size_t nodeIdx = 0; nodeIdx < nodeIndices.size(); ++nodeIdx)
    {
        if (nodeIndices[nodeIdx] >= allNodes_.size())
        {
            throw
                std::out_of_range
                (
                    "Node index "
                  + std::to_string(nodeIndices[nodeIdx])
                  + " out of range in skewness"
                  + " calculation for face "
                  + std::to_string(face.idx())
                );
        }

        Vector vertexToCentroid = 
            allNodes_[nodeIndices[nodeIdx]] - faceCentroid;

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
    const std::vector<size_t>& nodeIndices = face.nodeIndices();

    for (size_t i = 0; i < nodeIndices.size(); ++i)
    {
        if (nodeIndices[i] >= allNodes_.size())
        {
            throw
                std::out_of_range
                (
                    "Node index "
                  + std::to_string(nodeIndices[i])
                  + " out of range in boundary"
                  + " skewness calculation"
                  + " for face "
                  + std::to_string(face.idx())
                );
        }

        Vector vertexToCentroid = allNodes_[nodeIndices[i]] - faceCentroid;

        Scalar projection = std::abs(dot(skewnessDirection, vertexToCentroid));

        faceCharacteristicLength =
            std::max(faceCharacteristicLength, projection);
    }

    return skewnessMag / faceCharacteristicLength;
}


Scalar MeshChecker::calculateCellAspectRatio
(
    const Cell& cell
) const noexcept
{
    // Accumulate absolute face area components per direction
    Vector sumMagAreaComponents(0.0, 0.0, 0.0);

    const auto& faceIndices = cell.faceIndices();

    for (size_t i = 0; i < faceIndices.size(); ++i)
    {
        const Face& face = allFaces_[faceIndices[i]];
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
