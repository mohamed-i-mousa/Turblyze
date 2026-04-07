/******************************************************************************
 * @file Mesh.cpp
 * @brief Implementation of the owning Mesh container
 *****************************************************************************/

#include "Mesh.hpp"


Mesh::Mesh
(
    std::vector<Vector> nodes,
    std::vector<Face> faces,
    std::vector<Cell> cells,
    std::vector<BoundaryPatch> patches
)
    : nodes_(std::move(nodes)),
      faces_(std::move(faces)),
      cells_(std::move(cells)),
      patches_(std::move(patches))
{
    cellCount_ = cells_.size();
    faceCount_ = faces_.size();
}
