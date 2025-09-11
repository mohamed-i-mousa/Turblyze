/******************************************************************************
 * @file checkMesh.cpp
 * @brief Mesh quality assessment utilities
 * 
 * This file contains functions for checking mesh quality by analyzing
 * geometric properties such as face areas and cell volumes.
 *****************************************************************************/

#include "checkMesh.h"
#include "Scalar.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>

/**
 * @brief Calculate face non-orthogonality (OpenFOAM method)
 * @param ownCc Owner cell center
 * @param neiCc Neighbor cell center  
 * @param faceNormal Face normal vector (unit vector)
 * @return Cosine of angle between cell centers line and face normal
 */
static Scalar calculateFaceOrthogonality
(
    const Vector& ownerCellCentroid,
    const Vector& neighborCellCentroid,
    const Vector& faceNormal
)
{
    Vector d_PN = neighborCellCentroid - ownerCellCentroid;
    Scalar d_PN_mag = d_PN.magnitude();
    Scalar faceNormal_mag = faceNormal.magnitude();

    Scalar cosAngle = dot(d_PN, faceNormal) / (d_PN_mag * faceNormal_mag + vSmallValue);
    return std::clamp(cosAngle, S(-1.0), S(1.0));
}

/**
 * @brief Calculate face skewness (OpenFOAM intersection point method)
 * @param ownCc Owner cell center
 * @param neiCc Neighbor cell center
 * @param faceCtr Face center
 * @param faceNormal Face normal vector
 * @return Skewness value (0 = perfect, >1 = skewed)
 */
static Scalar calculateFaceSkewness
(
    const Vector& ownerCellCentroid,
    const Vector& neighborCellCentroid, 
    const Vector& faceCentroid,
    const Vector& faceNormal
)
{
    
    Vector d_Pf = faceCentroid - ownerCellCentroid;

    Vector d_PN = neighborCellCentroid - ownerCellCentroid;
    
    // Skewness vector (OpenFOAM method)
    Vector skewnessVector = d_Pf - ((dot(faceNormal, d_Pf)) / (dot(faceNormal, d_PN) + vSmallValue)) * d_PN;
    
    // Normalization distance (simplified without face vertices)
    Scalar normalizationDistance = S(0.2) * d_PN.magnitude() + vSmallValue;
    
    // Return normalized skewness
    return skewnessVector.magnitude() / normalizationDistance; 
}

/**
 * @brief Calculate boundary face skewness  
 * @param cellCc Cell center
 * @param faceCtr Face center
 * @param faceNormal Face normal vector
 * @return Boundary skewness value
 */
static Scalar calculateBoundarySkewness
(
    const Vector& ownerCellCentroid,
    const Vector& faceCentroid,
    const Vector& faceNormal
)
{
    // Cell-to-face vector (similar to internal face Cpf)
    Vector d_Pf = faceCentroid - ownerCellCentroid;
    
    // Skewness vector: tangential component (deviation from normal line)
    Vector skewnessVector = d_Pf - dot(faceNormal, d_Pf) * faceNormal;
    
    // Normalization distance 
    Scalar normalizationDistance = d_Pf.magnitude();
    if (normalizationDistance < vSmallValue)
    {
        return S(0.0);
    }
    
    // Return normalized skewness magnitude
    return skewnessVector.magnitude() / normalizationDistance;
}

/**
 * @brief Calculate cell aspect ratio (OpenFOAM directional method)
 * @param cell Cell object
 * @param allFaces Vector of all faces
 * @return Aspect ratio (1.0 = perfect cube, higher = elongated)
 */
static Scalar calculateCellAspectRatio
(
    const Cell& cell,
    const std::vector<Face>& allFaces
)
{
    // Accumulate face area components in each direction
    Vector sumMagAreaComponents(0.0, 0.0, 0.0);
    
    const auto& faceIndices = cell.faceIndices();
    const auto& faceSigns = cell.faceSigns();
    
    for (size_t i = 0; i < faceIndices.size(); ++i)
    {
        const Face& face = allFaces[faceIndices[i]];
        Vector areaVec = face.normal() * face.area();
        
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
        ({
            sumMagAreaComponents.x(), 
            sumMagAreaComponents.y(), 
            sumMagAreaComponents.z()
        });

    Scalar maxComponent = 
        std::max
        ({
            sumMagAreaComponents.x(),
            sumMagAreaComponents.y(),
            sumMagAreaComponents.z()
        });
    
    if (minComponent < vSmallValue)
    {
        return S(1000.0); // Very high aspect ratio for degenerate cells
    }
    
    Scalar directionalAspect = maxComponent / minComponent;
    
    // Add hydraulic aspect ratio for 3D cells (OpenFOAM method)
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

void checkMesh
(
    const std::vector<Face>& allFaces,
    const std::vector<Cell>& allCells
)
{
    std::cout << "\n--- Mesh Quality Check ---" << std::endl;
    
    // Initialize with first face/cell values
    if (allFaces.empty() || allCells.empty())
    {
        std::cout << "Warning: Empty mesh detected!" << std::endl;
        return;
    }
    
    // Face area statistics
    Scalar minFaceArea = allFaces[0].area();
    Scalar maxFaceArea = allFaces[0].area();
    size_t minFaceId = allFaces[0].id();
    size_t maxFaceId = allFaces[0].id();
    
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
    
    // Process all faces
    for (const auto& face : allFaces)
    {
        // Area statistics
        if (face.area() < minFaceArea)
        {
            minFaceArea = face.area();
            minFaceId = face.id();
        }
        if (face.area() > maxFaceArea)
        {
            maxFaceArea = face.area();
            maxFaceId = face.id();
        }
        
        // Check for faces with area smaller than 1e-12
        if (face.area() < minArea)
        {
            smallAreaFaces.push_back(face.id());
        }
        
        // Calculate non-orthogonality and skewness
        if (face.isBoundary())
        {
            // Boundary face: calculate skewness only
            const Cell& ownerCell = allCells[face.ownerCell()];
            Scalar skew = 
                calculateBoundarySkewness
                (
                    ownerCell.centroid(),
                    face.centroid(),
                    face.normal()
                );
            
            if (skew > maxSkewness)
            {
                maxSkewness = skew;
                maxSkewFaceId = face.id();
            }
            if (skew > S(4.0))
            {
                highSkewFaces.push_back(face.id());
            }
        }
        else
        {
            // Internal face: calculate both non-orthogonality and skewness
            const Cell& ownerCell = allCells[face.ownerCell()];
            const Cell& neighborCell = allCells[face.neighborCell().value()];
            
            // Non-orthogonality (angle in degrees)
            Scalar ortho = 
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
                maxNonOrthFaceId = face.id();
            }
            
            // Threshold: 70 degrees
            if (angleDeg > S(70.0))
            {
                severeNonOrthFaces.push_back(face.id());
            }
            
            // Skewness
            Scalar skew = calculateFaceSkewness(
                ownerCell.centroid(),
                neighborCell.centroid(),
                face.centroid(),
                face.normal()
            );
            
            if (skew > maxSkewness)
            {
                maxSkewness = skew;
                maxSkewFaceId = face.id();
            }
            
            // OpenFOAM threshold: 4.0
            if (skew > S(4.0))
            {
                highSkewFaces.push_back(face.id());
            }
        }
    }
    
    // Calculate average non-orthogonality
    Scalar avgNonOrthogonality = 
        std::acos(totalNonOrthogonality / S(nonOrthCount));
    
    avgNonOrthogonality = avgNonOrthogonality * S(180.0) / S(M_PI);

    // Cell volume and aspect ratio statistics
    Scalar minCellVolume = allCells[0].volume();
    Scalar maxCellVolume = allCells[0].volume();
    size_t minCellId = allCells[0].id();
    size_t maxCellId = allCells[0].id();
    
    Scalar maxAspectRatio = 0.0;
    size_t maxAspectCellId = 0;
    std::vector<size_t> highAspectCells; // > 100
    
    std::vector<size_t> smallVolumeCells;
    
    for (const auto& cell : allCells)
    {
        // Volume statistics
        if (cell.volume() < minCellVolume)
        {
            minCellVolume = cell.volume();
            minCellId = cell.id();
        }
        if (cell.volume() > maxCellVolume)
        {
            maxCellVolume = cell.volume();
            maxCellId = cell.id();
        }
        
        // Check for cells with volume smaller than 1e-30
        if (cell.volume() < minVolume)
        {
            smallVolumeCells.push_back(cell.id());
        }
        
        // Calculate aspect ratio
        Scalar aspectRatio = calculateCellAspectRatio(cell, allFaces);
        if (aspectRatio > maxAspectRatio)
        {
            maxAspectRatio = aspectRatio;
            maxAspectCellId = cell.id();
        }
        
        // High aspect ratio threshold
        if (aspectRatio > S(100.0))
        {
            highAspectCells.push_back(cell.id());
        }
    }
    
    // Store current precision setting
    std::streamsize oldPrecision = std::cout.precision();
    
    // Report results with scientific notation for small values
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
    
    // Restore original precision settings
    std::cout << std::setprecision(oldPrecision);
    
    // Quality warnings for small areas/volumes
    if (!smallAreaFaces.empty())
    {
        std::cout   << "\nQuality Check - Small Face Areas:" << std::endl;
        std::cout   << "  Found " << smallAreaFaces.size() 
                    << " faces with area < "
                    << std::scientific << std::setprecision(0) << minArea 
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
                    << minVolume << " m³" << std::endl;
        
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