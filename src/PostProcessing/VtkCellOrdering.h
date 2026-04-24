/******************************************************************************
 * @file VtkCellOrdering.h
 * @brief Internal helpers for reconstructing VTK cell node orderings
 *
 * @details Topology utilities that translate face-based cell connectivity
 * (as stored in the internal mesh) into the per-vertex orderings that VTK
 * expects for hexahedra, wedges, and pyramids. Also provides the small
 * geometric helpers (face centroid and normal) used to disambiguate
 * opposite faces and winding directions.
 *
 * @note These declarations are intentionally not exported through the
 * public VtkWriter API — they are consumed only by VtkWriter.cpp.
 *****************************************************************************/

#pragma once

#include <span>
#include <unordered_set>
#include <vector>

#include <vtkType.h>

#include "Vector.h"


namespace VTK
{

/**
 * @brief Compute the arithmetic centroid of a polygonal face
 * @param faceNodes Indices into allNodes for this face
 * @param allNodes Mesh node coordinates
 * @return Face centroid
 */
[[nodiscard]] Vector faceCentroid
(
    std::span<const size_t> faceNodes,
    std::span<const Vector> allNodes
);

/**
 * @brief Compute an unnormalized face normal from the first three vertices
 * @param faceNodes Indices into allNodes for this face
 * @param allNodes Mesh node coordinates
 * @return Cross product of edges (v1-v0) x (v2-v0)
 */
[[nodiscard]] Vector faceNormal
(
    std::span<const size_t> faceNodes,
    std::span<const Vector> allNodes
);

/**
 * @brief Produce VTK_HEXAHEDRON node ordering from face connectivity
 * @param faceNodeLists Per-face lists of node indices
 * @param uniqueNodes All distinct node indices in the cell
 * @param allNodes Mesh node coordinates
 * @return 8-element VTK-ordered node list, or empty on failure
 */
[[nodiscard]] std::vector<vtkIdType> orderHexahedronNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::unordered_set<size_t>& uniqueNodes,
    std::span<const Vector> allNodes
);

/**
 * @brief Produce VTK_WEDGE node ordering from face connectivity
 * @param faceNodeLists Per-face lists of node indices
 * @param uniqueNodes All distinct node indices in the cell
 * @param allNodes Mesh node coordinates
 * @return 6-element VTK-ordered node list, or empty on failure
 */
[[nodiscard]] std::vector<vtkIdType> orderWedgeNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::unordered_set<size_t>& uniqueNodes,
    std::span<const Vector> allNodes
);

/**
 * @brief Produce VTK_PYRAMID node ordering from face connectivity
 * @param faceNodeLists Per-face lists of node indices
 * @param uniqueNodes All distinct node indices in the cell
 * @param allNodes Mesh node coordinates
 * @return 5-element VTK-ordered node list, or empty on failure
 */
[[nodiscard]] std::vector<vtkIdType> orderPyramidNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::unordered_set<size_t>& uniqueNodes,
    std::span<const Vector> allNodes
);

} // namespace VTK
