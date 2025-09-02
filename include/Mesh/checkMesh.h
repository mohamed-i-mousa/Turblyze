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