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
#include <iostream>
#include <iomanip>
#include <algorithm>

void performMeshCheck(const std::vector<Face>& allFaces, const std::vector<Cell>& allCells)
{
    std::cout << "\n--- Mesh Quality Check ---" << std::endl;
    
    // Initialize with first face/cell values
    if (allFaces.empty() || allCells.empty()) {
        std::cout << "Warning: Empty mesh detected!" << std::endl;
        return;
    }
    
    // Face area statistics
    Scalar minFaceArea = allFaces[0].area;
    Scalar maxFaceArea = allFaces[0].area;
    size_t minFaceId = allFaces[0].id;
    size_t maxFaceId = allFaces[0].id;
    
    for (const auto& face : allFaces) {
        if (face.area < minFaceArea) {
            minFaceArea = face.area;
            minFaceId = face.id;
        }
        if (face.area > maxFaceArea) {
            maxFaceArea = face.area;
            maxFaceId = face.id;
        }
    }
    
    // Cell volume statistics  
    Scalar minCellVolume = allCells[0].volume;
    Scalar maxCellVolume = allCells[0].volume;
    size_t minCellId = allCells[0].id;
    size_t maxCellId = allCells[0].id;
    
    for (const auto& cell : allCells) {
        if (cell.volume < minCellVolume) {
            minCellVolume = cell.volume;
            minCellId = cell.id;
        }
        if (cell.volume > maxCellVolume) {
            maxCellVolume = cell.volume;
            maxCellId = cell.id;
        }
    }
    
    // Calculate aspect ratios
    Scalar areaRatio = maxFaceArea / minFaceArea;
    Scalar volumeRatio = maxCellVolume / minCellVolume;
    
    // Store current precision setting
    std::streamsize oldPrecision = std::cout.precision();
    
    // Report results with scientific notation for small values
    std::cout << "Face Area Statistics:" << std::endl;
    std::cout << "  Minimum area: " << std::scientific << std::setprecision(6) 
              << minFaceArea << " m² (face " << minFaceId << ")" << std::endl;
    std::cout << "  Maximum area: " << std::scientific << std::setprecision(6) 
              << maxFaceArea << " m² (face " << maxFaceId << ")" << std::endl;
    std::cout << "  Area ratio (max/min): " << std::fixed << std::setprecision(2) 
              << areaRatio << std::endl;
    
    std::cout << "Cell Volume Statistics:" << std::endl;
    std::cout << "  Minimum volume: " << std::scientific << std::setprecision(6) 
              << minCellVolume << " m³ (cell " << minCellId << ")" << std::endl;
    std::cout << "  Maximum volume: " << std::scientific << std::setprecision(6) 
              << maxCellVolume << " m³ (cell " << maxCellId << ")" << std::endl;
    std::cout << "  Volume ratio (max/min): " << std::fixed << std::setprecision(2) 
              << volumeRatio << std::endl;
    
    // Restore original precision settings
    std::cout << std::setprecision(oldPrecision);
    
    // Quality warnings
    if (areaRatio > 1e6) {
        std::cout << "Warning: Large face area variation may cause numerical issues!" << std::endl;
    }
    if (volumeRatio > 1e6) {
        std::cout << "Warning: Large cell volume variation may cause numerical issues!" << std::endl;
    }
    if (minFaceArea < 1e-12) {
        std::cout << "Warning: Very small face areas detected!" << std::endl;
    }
    if (minCellVolume < 1e-18) {
        std::cout << "Warning: Very small cell volumes detected!" << std::endl;
    }
}