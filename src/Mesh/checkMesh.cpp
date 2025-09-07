/**
 * @file checkMesh.cpp
 * @brief Mesh quality assessment utilities
 * 
 * This file contains functions for checking mesh quality by analyzing
 * geometric properties such as face areas and cell volumes.
 * 
 * @author Mohamed Mousa
 * @date 2025
 */

#include "checkMesh.h"
#include "Scalar.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

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
    
    for (const auto& face : allFaces)
    {
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
    }
    
    // Cell volume statistics  
    Scalar minCellVolume = allCells[0].volume();
    Scalar maxCellVolume = allCells[0].volume();
    size_t minCellId = allCells[0].id();
    size_t maxCellId = allCells[0].id();
    
    std::vector<size_t> smallVolumeCells;
    
    for (const auto& cell : allCells)
    {
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
    }
    
    // Store current precision setting
    std::streamsize oldPrecision = std::cout.precision();
    
    // Report results with scientific notation for small values
    std::cout << "Face Area Statistics:" << std::endl;

    std::cout   << "  Minimum area: " << std::scientific 
                << std::setprecision(6) << minFaceArea << " m² (face "
                << minFaceId << ")" << std::endl;

    std::cout   << "  Maximum area: " << std::scientific
                << std::setprecision(6) << maxFaceArea << " m² (face " 
                << maxFaceId << ")" << std::endl;
    
    std::cout   << "Cell Volume Statistics:" << std::endl;

    std::cout   << "  Minimum volume: " << std::scientific
                << std::setprecision(6) << minCellVolume << " m³ (cell " 
                << minCellId << ")" << std::endl;

    std::cout   << "  Maximum volume: " << std::scientific 
                << std::setprecision(6) << maxCellVolume << " m³ (cell "
                << maxCellId << ")" << std::endl;
    
    // Restore original precision settings
    std::cout << std::setprecision(oldPrecision);
    
    // Quality warnings and detailed reporting
    if (!smallAreaFaces.empty())
    {
        std::cout << "\nQuality Check - Small Face Areas:" << std::endl;
        std::cout << "  Found " << smallAreaFaces.size() << " faces with area < " 
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
        std::cout << "\nQuality Check - Small Cell Volumes:" << std::endl;
        std::cout << "  Found " << smallVolumeCells.size() << " cells with volume < " 
                  << std::scientific << std::setprecision(0) << minVolume 
                  << " m³" << std::endl;
        
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
    
    if (smallAreaFaces.empty() && smallVolumeCells.empty())
    {
        std::cout << "\n✓ All faces and cells pass minimum size requirements" << std::endl;
    }
}