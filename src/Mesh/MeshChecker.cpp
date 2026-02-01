/******************************************************************************
 * @file checkMesh.cpp
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
)
:
    allNodes_(nodes),
    allFaces_(faces),
    allCells_(cells)
{}


// ***************************** Private Methods ******************************

Scalar MeshChecker::calculateFaceOrthogonality
(
    const Vector& ownerCellCentroid,
    const Vector& neighborCellCentroid,
    const Vector& faceNormal
) const
{
    Vector d_PN = neighborCellCentroid - ownerCellCentroid;
    Scalar d_PN_mag = d_PN.magnitude();
    Scalar faceNormal_mag = faceNormal.magnitude();

    Scalar cosAngle =
        dot(d_PN, faceNormal) / (d_PN_mag * faceNormal_mag + vSmallValue);

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
    Vector d_Pf = faceCentroid - ownerCellCentroid;
    Vector d_PN = neighborCellCentroid - ownerCellCentroid;

    Vector skewnessVector =
        d_Pf
      - ((dot(faceNormal, d_Pf)) / (dot(faceNormal, d_PN) + vSmallValue))
      * d_PN;

    Scalar skewnessMag = skewnessVector.magnitude();

    Vector skewnessDirection = skewnessVector / (skewnessMag + smallValue);

    // Characteristic face dimension: start with empirical approximation
    Scalar faceCharacteristicLength = S(0.2) * d_PN.magnitude() + vSmallValue;

    // Refine by finding maximum vertex extent in skewness direction
    const std::vector<size_t>& nodeIndices = face.nodeIndices();
    for (size_t nodeIdx = 0; nodeIdx < nodeIndices.size(); ++nodeIdx)
    {
        if (nodeIndices[nodeIdx] >= allNodes_.size())
        {
            throw   std::out_of_range
                    (
                        "Node index " + std::to_string(nodeIndices[nodeIdx])
                      + " out of range in skewness calculation for face "
                      + std::to_string(face.idx())
                    );
        }

        Vector vertexToCentroid = allNodes_[nodeIndices[nodeIdx]] - faceCentroid;

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
    Vector d_Pf = faceCentroid - ownerCellCentroid;

    // Virtual d_PN for boundary faces
    Vector d_PN = S(2.0) * dot(faceNormal, d_Pf) * faceNormal;

    Vector skewnessVector = d_Pf - dot(faceNormal, d_Pf) * faceNormal;

    Scalar skewnessMag = skewnessVector.magnitude();

    Vector skewnessDirection = skewnessVector / (skewnessMag + smallValue);

    // Characteristic face dimension: start with empirical approximation
    Scalar faceCharacteristicLength = S(0.4) * d_PN.magnitude() + vSmallValue;

    // Refine by finding maximum vertex extent in skewness direction
    const std::vector<size_t>& nodeIndices = face.nodeIndices();
    for (size_t i = 0; i < nodeIndices.size(); ++i)
    {
        if (nodeIndices[i] >= allNodes_.size())
        {
            throw   std::out_of_range
                    (
                        "Node index " + std::to_string(nodeIndices[i])
                      + " out of range in boundary skewness calculation "
                      + "calculation for face "
                      + std::to_string(face.idx())
                    );
        }

        Vector vertexToCentroid = allNodes_[nodeIndices[i]] - faceCentroid;

        Scalar projection = 
            std::abs(dot(skewnessDirection, vertexToCentroid));

        faceCharacteristicLength = 
            std::max(faceCharacteristicLength, projection);
    }

    return skewnessMag / faceCharacteristicLength;
}


Scalar MeshChecker::calculateCellAspectRatio(const Cell& cell) const
{
    // Accumulate face area components in each direction
    Vector sumMagAreaComponents(0.0, 0.0, 0.0);

    const auto& faceIndices = cell.faceIndices();
    const auto& faceSigns = cell.faceSigns();

    for (size_t i = 0; i < faceIndices.size(); ++i)
    {
        const Face& face = allFaces_[faceIndices[i]];
        Vector areaVec = face.normal() * face.projectedArea();

        // Account for face orientation relative to cell
        if (faceSigns[i] < 0)
        {
            areaVec = areaVec * S(-1.0);
        }

        // Accumulate absolute components
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
    Scalar minComponent =
                std::min
                (
                    {
                        sumMagAreaComponents.x(),
                        sumMagAreaComponents.y(),
                        sumMagAreaComponents.z()
                    }
                );

    Scalar  maxComponent =
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
    Scalar totalSurfaceArea =
        sumMagAreaComponents.x()
      + sumMagAreaComponents.y()
      + sumMagAreaComponents.z();

    Scalar volume = cell.volume();

    if (volume > vSmallValue)
    {
        // Hydraulic aspect ratio: (1/6) * A / V^(2/3)
        Scalar hydraulicAspect =
            (S(1.0)/S(6.0)) * totalSurfaceArea
          / std::pow(volume, S(2.0)/S(3.0));

        directionalAspect = std::max(directionalAspect, hydraulicAspect);
    }

    return directionalAspect;
}


// ****************************** Public Methods ******************************

void MeshChecker::check() const
{
    std::cout << "\n--- Mesh Quality Check ---" << std::endl;

    if (allFaces_.empty() || allCells_.empty())
    {
        std::cout << "Warning: Empty mesh detected!" << std::endl;
        return;
    }

    // Face area statistics
    Scalar minFaceArea = allFaces_[0].projectedArea();
    Scalar maxFaceArea = allFaces_[0].projectedArea();
    size_t minFaceId = allFaces_[0].idx();
    size_t maxFaceId = allFaces_[0].idx();

    std::vector<size_t> smallAreaFaces;

    // Non-orthogonality statistics
    Scalar maxNonOrthogonality = 0.0;
    Scalar totalNonOrthogonality = 0.0;
    size_t nonOrthCount = 0;
    size_t maxNonOrthFaceId = 0;
    std::vector<size_t> severeNonOrthFaces; // > 70 degrees

    // Skewness statistics
    Scalar maxSkewness = 0.0;
    size_t maxSkewFaceId = 0;
    std::vector<size_t> highSkewFaces; // > 4.0

    for (const auto& face : allFaces_)
    {
        // Area statistics
        if (face.projectedArea() < minFaceArea)
        {
            minFaceArea = face.projectedArea();
            minFaceId = face.idx();
        }
        if (face.projectedArea() > maxFaceArea)
        {
            maxFaceArea = face.projectedArea();
            maxFaceId = face.idx();
        }

        if (face.projectedArea() < minArea_)
        {
            smallAreaFaces.push_back(face.idx());
        }

        // Calculate non-orthogonality and skewness
        if (face.isBoundary())
        {
            // Boundary face: calculate skewness only
            if (face.ownerCell() >= allCells_.size())
            {
                throw   std::out_of_range
                        (
                            "Owner cell index " 
                          + std::to_string(face.ownerCell())
                          + " out of range for face " 
                          + std::to_string(face.idx())
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
                maxSkewFaceId = face.idx();
            }
            if (skew > S(4.0))
            {
                highSkewFaces.push_back(face.idx());
            }
        }
        else
        {
            // Internal face: calculate both non-orthogonality and skewness
            if (face.ownerCell() >= allCells_.size())
            {
                throw   std::out_of_range
                        (
                            "Owner cell index " 
                          + std::to_string(face.ownerCell())
                          + " out of range for face " 
                          + std::to_string(face.idx())
                        );
            }

            if (face.neighborCell().value() >= allCells_.size())
            {
                throw   std::out_of_range
                        (
                            "Neighbor cell index "
                          + std::to_string(face.neighborCell().value())
                          + " out of range for face " 
                          + std::to_string(face.idx())
                        );
            }

            const Cell& ownerCell = allCells_[face.ownerCell()];
            const Cell& neighborCell = allCells_[face.neighborCell().value()];

            // Non-orthogonality (angle in degrees)
            Scalar  ortho =
                        calculateFaceOrthogonality
                        (
                            ownerCell.centroid(),
                            neighborCell.centroid(),
                            face.normal()
                        );

            // Convert to angle in degrees (with safety check)
            Scalar angleRad = std::acos(ortho);
            Scalar angleDeg = angleRad * S(180.0) / S(M_PI);

            totalNonOrthogonality += ortho;

            nonOrthCount++;

            if (angleDeg > maxNonOrthogonality)
            {
                maxNonOrthogonality = angleDeg;
                maxNonOrthFaceId = face.idx();
            }

            if (angleDeg > S(70.0))
            {
                severeNonOrthFaces.push_back(face.idx());
            }

            // Skewness
            Scalar  skew =   
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
                maxSkewFaceId = face.idx();
            }

            if (skew > S(4.0))
            {
                highSkewFaces.push_back(face.idx());
            }
        }
    }

    // Calculate average non-orthogonality
    Scalar avgNonOrthogonality =
        std::acos(totalNonOrthogonality / S(nonOrthCount));

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
        if (aspectRatio > S(100.0))
        {
            highAspectCells.push_back(cell.idx());
        }
    }

    // Store current format flags and precision
    std::ios_base::fmtflags oldFlags = std::cout.flags();
    std::streamsize oldPrecision = std::cout.precision();

    std::cout << "\nFace Area Statistics:" << std::endl;

    std::cout   << "  Minimum area: " << std::scientific
                << std::setprecision(6) << minFaceArea << " m² (face "
                << minFaceId << ")" << std::endl;

    std::cout   << "  Maximum area: " << std::scientific
                << std::setprecision(6) << maxFaceArea << " m² (face "
                << maxFaceId << ")" << std::endl;

    std::cout   << "\nCell Volume Statistics:" << std::endl;

    std::cout   << "  Minimum volume: " << std::scientific
                << std::setprecision(6) << minCellVolume << " m³ (cell "
                << minCellId << ")" << std::endl;

    std::cout   << "  Maximum volume: " << std::scientific
                << std::setprecision(6) << maxCellVolume << " m³ (cell "
                << maxCellId << ")" << std::endl;

    // Non-orthogonality statistics
    std::cout << "\nNon-Orthogonality Statistics:" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Maximum: " << maxNonOrthogonality << "° (face "
              << maxNonOrthFaceId << ")" << std::endl;
    std::cout << "  Average: " << avgNonOrthogonality << "°" << std::endl;

    if (!severeNonOrthFaces.empty())
    {
        std::cout << "  WARNING: " << severeNonOrthFaces.size()
                  << " faces with non-orthogonality > 70°" << std::endl;
        if (severeNonOrthFaces.size() <= 10)
        {
            std::cout << "  Face IDs: ";
            for (size_t i = 0; i < severeNonOrthFaces.size(); ++i)
            {
                std::cout << severeNonOrthFaces[i];
                if (i < severeNonOrthFaces.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
    }

    // Skewness statistics
    std::cout << "\nSkewness Statistics:" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  Maximum: " << maxSkewness << " (face "
              << maxSkewFaceId << ")" << std::endl;

    if (!highSkewFaces.empty())
    {
        std::cout << "  WARNING: " << highSkewFaces.size()
                  << " faces with skewness > 4.0" << std::endl;
        if (highSkewFaces.size() <= 10)
        {
            std::cout << "  Face IDs: ";
            for (size_t i = 0; i < highSkewFaces.size(); ++i)
            {
                std::cout << highSkewFaces[i];
                if (i < highSkewFaces.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
    }

    // Aspect ratio statistics
    std::cout << "\nAspect Ratio Statistics:" << std::endl;
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  Maximum: " << maxAspectRatio << " (cell "
              << maxAspectCellId << ")" << std::endl;

    if (!highAspectCells.empty())
    {
        std::cout << "  WARNING: " << highAspectCells.size()
                  << " cells with aspect ratio > 100" << std::endl;
        if (highAspectCells.size() <= 10)
        {
            std::cout << "  Cell IDs: ";
            for (size_t i = 0; i < highAspectCells.size(); ++i)
            {
                std::cout << highAspectCells[i];
                if (i < highAspectCells.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
    }

    // Restore original format flags and precision
    std::cout.flags(oldFlags);
    std::cout.precision(oldPrecision);

    // Quality warnings for small areas/volumes
    if (!smallAreaFaces.empty())
    {
        std::cout   << "\nQuality Check - Small Face Areas:" << std::endl;
        std::cout   << "  Found " << smallAreaFaces.size()
                    << " faces with area < "
                    << std::scientific << std::setprecision(0) << minArea_
                    << " m²" << std::endl;

        if (smallAreaFaces.size() <= 10)
        {
            std::cout << "  Face IDs: ";
            for (size_t i = 0; i < smallAreaFaces.size(); ++i)
            {
                std::cout << smallAreaFaces[i];
                if (i < smallAreaFaces.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
        else
        {
            std::cout << "  First 10 face IDs: ";
            for (size_t i = 0; i < 10; ++i)
            {
                std::cout << smallAreaFaces[i];
                if (i < 9) std::cout << ", ";
            }
            std::cout << " ..." << std::endl;
        }
    }

    if (!smallVolumeCells.empty())
    {
        std::cout   << "\nQuality Check - Small Cell Volumes:" << std::endl;
        std::cout   << "  Found " << smallVolumeCells.size()
                    << " cells with volume < "
                    << std::scientific << std::setprecision(0)
                    << minVolume_ << " m³" << std::endl;

        if (smallVolumeCells.size() <= 10)
        {
            std::cout << "  Cell IDs: ";
            for (size_t i = 0; i < smallVolumeCells.size(); ++i)
            {
                std::cout << smallVolumeCells[i];
                if (i < smallVolumeCells.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
        else
        {
            std::cout << "  First 10 cell IDs: ";
            for (size_t i = 0; i < 10; ++i)
            {
                std::cout << smallVolumeCells[i];
                if (i < 9) std::cout << ", ";
            }
            std::cout << " ..." << std::endl;
        }
    }

    // Overall mesh quality summary
    std::cout << "\n--- Mesh Quality Summary ---" << std::endl;

    bool goodQuality = true;

    if (maxNonOrthogonality > S(70.0))
    {
        std::cout << "⚠ Non-orthogonality exceeds 70° threshold" << std::endl;
        goodQuality = false;
    }

    if (maxSkewness > S(4.0))
    {
        std::cout << "⚠ Skewness exceeds 4.0 threshold" << std::endl;
        goodQuality = false;
    }

    if (maxAspectRatio > S(100.0))
    {
        std::cout << "⚠ Aspect ratio exceeds 100 threshold" << std::endl;
        goodQuality = false;
    }

    if (smallAreaFaces.empty() && smallVolumeCells.empty() && goodQuality)
    {
        std::cout   << "✓ All mesh quality metrics within acceptable ranges"
                    << std::endl;
    }
}