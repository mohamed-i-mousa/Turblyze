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
 * - Geometric property validation for numerical stability
 * - Scientific notation output for proper display of small values
 * 
 * The mesh checker helps identify:
 * - Degenerate faces or cells that could cause solver instability
 * - Extreme aspect ratios that might require special treatment
 * - Volume conservation issues and geometric inconsistencies
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
 * Reports minimum and maximum face areas and cell volumes to help
 * assess mesh quality and identify potential numerical issues.
 * Uses scientific notation to properly display very small values.
 */
void checkMesh
(
    const std::vector<Face>& allFaces,
    const std::vector<Cell>& allCells
);

#endif