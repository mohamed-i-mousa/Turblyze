/******************************************************************************
 * @file Mesh.hpp
 * @brief Owning mesh container and lightweight view provider
 *
 * @details This header defines the Mesh class, which owns all mesh data
 * (nodes, faces, cells, boundary patches) and provides span-based views
 * to consumers. It also manages cell and face counts used by field
 * containers at construction time.
 *
 * @class Mesh
 * - Owns nodes, faces, cells, and boundary patches via std::vector members
 * - Const span accessors for read-only consumers (const Mesh&)
 * - Mutable span accessors for the mesh preparation phase
 * - Static registration and retrieval of cell and face counts
 *****************************************************************************/

#pragma once

#include <span>
#include <cstddef>
#include <vector>

#include "Vector.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "BoundaryPatch.hpp"


class Mesh
{
public:

    /// Construct empty mesh (for delayed initialization)
    Mesh() = default;

    /**
     * @brief Construct mesh by taking ownership of data vectors
     * @param nodes  Node coordinate data (moved in)
     * @param faces  Face topology data (moved in)
     * @param cells  Cell topology data (moved in)
     * @param patches Boundary patch data (moved in)
     */
    Mesh
    (
        std::vector<Vector> nodes,
        std::vector<Face> faces,
        std::vector<Cell> cells,
        std::vector<BoundaryPatch> patches
    );

    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    Mesh(Mesh&&) = default;
    Mesh& operator=(Mesh&&) = default;

// Const accessor methods (for const Mesh& consumers)

    /// Read-only view of node coordinates
    std::span<const Vector> nodes() const noexcept
    {
        return nodes_;
    }

    /// Read-only view of faces
    std::span<const Face> faces() const noexcept
    {
        return faces_;
    }

    /// Read-only view of cells
    std::span<const Cell> cells() const noexcept
    {
        return cells_;
    }

    /// Read-only view of boundary patches
    std::span<const BoundaryPatch> patches() const noexcept
    {
        return patches_;
    }

// Mutable accessor methods (for mesh preparation phase only)

    /// Mutable view of faces — use only during prepareMesh()
    std::span<Face> faces() noexcept { return faces_; }

    /// Mutable view of cells — use only during prepareMesh()
    std::span<Cell> cells() noexcept { return cells_; }

// Size accessor methods

    /// Number of nodes in the mesh
    size_t numNodes() const noexcept { return nodes_.size(); }

    /// Number of faces in the mesh
    size_t numFaces() const noexcept { return faces_.size(); }

    /// Number of cells in the mesh
    size_t numCells() const noexcept { return cells_.size(); }

    /// Cell count registered at startup (used by CellData/FaceData)
    static size_t cellCount() noexcept { return cellCount_; }

    /// Face count registered at startup (used by CellData/FaceData)
    static size_t faceCount() noexcept { return faceCount_; }

private:

    /// Mesh node coordinates
    std::vector<Vector> nodes_;

    /// Mesh face topology
    std::vector<Face> faces_;

    /// Mesh cell topology
    std::vector<Cell> cells_;

    /// Boundary patch descriptors
    std::vector<BoundaryPatch> patches_;

    /// Cell count for field container initialization (set once at startup)
    static inline size_t cellCount_ = 0;

    /// Face count for field container initialization (set once at startup)
    static inline size_t faceCount_ = 0;
};
