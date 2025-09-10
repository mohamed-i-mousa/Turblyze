/******************************************************************************
 * @file checkMesh.h
 * @brief Mesh quality assessment and diagnostic utilities
 * 
 * This header provides mesh quality checking functions that analyze geometric
 * properties and report statistics to help identify potential numerical 
 * issues.
 * 
 * The diagnostics include face area distribution, cell volume analysis, and
 * connectivity validation to ensure mesh suitability for CFD calculations.
 * 
 * Key mesh quality metrics:
 * - Minimum and maximum face areas with statistics
 * - Cell volume distribution and extrema detection
 * 
 *****************************************************************************/

#ifndef CHECKMESH_H
#define CHECKMESH_H

#include <vector>
#include "Face.h"
#include "Cell.h"

/**
 * @brief Performs mesh quality checks and reports statistics
 * @param allFaces Vector containing all mesh faces
 * @param allCells Vector containing all mesh cells
 * 
 * Reports mesh quality statistics including:
 * - Face areas and cell volumes (min/max with scientific notation)
 */
void checkMesh
(
    const std::vector<Face>& allFaces,
    const std::vector<Cell>& allCells
);

#endif