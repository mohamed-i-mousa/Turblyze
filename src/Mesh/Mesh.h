/******************************************************************************
 * @file Mesh.h
 * @brief Owning mesh data and lightweight view provider
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
 * - Static retrieval of cell and face counts
 *****************************************************************************/

#pragma once

#include <span>
#include <cstddef>
#include <vector>

#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"


class Mesh
{
public:

    /// Construct empty mesh
    Mesh() = default;

    /**
     * @brief Construct mesh by taking ownership of data vectors
     * @param nodes  Node coordinate data
     * @param faces  Face topology data
     * @param cells  Cell topology data
     * @param patches Boundary patch data
     */
    Mesh
    (
        std::vector<Vector> nodes,
        std::vector<Face> faces,
        std::vector<Cell> cells,
        std::vector<BoundaryPatch> patches
    )
    : 
        nodes_(std::move(nodes)),
        faces_(std::move(faces)),
        cells_(std::move(cells)),
        patches_(std::move(patches))
    {
        cellCount_ = cells_.size();
        faceCount_ = faces_.size();
    }

    /// Copy constructor and assignment - Not copyable (copy is expensive)
    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    /// Move constructor and assignment
    Mesh(Mesh&&) noexcept = default;
    Mesh& operator=(Mesh&&) noexcept = default;

    /// Destructor
    ~Mesh() noexcept = default;

// Const accessor methods

    /// Read-only view of node coordinates
    [[nodiscard]] std::span<const Vector> nodes() const noexcept
    {
        return nodes_;
    }

    /// Read-only view of faces
    [[nodiscard]] std::span<const Face> faces() const noexcept
    {
        return faces_;
    }

    /// Read-only view of cells
    [[nodiscard]] std::span<const Cell> cells() const noexcept
    {
        return cells_;
    }

    /// Read-only view of boundary patches
    [[nodiscard]] std::span<const BoundaryPatch> patches() const noexcept
    {
        return patches_;
    }

// Mutable accessor methods (during prepareMesh() only)

    /// Mutable view of faces
    [[nodiscard]] std::span<Face> faces() noexcept
    {
        return faces_;
    }

    /// Mutable view of cells
    [[nodiscard]] std::span<Cell> cells() noexcept
    {
        return cells_;
    }

// Size accessor methods

    /// Number of nodes in the mesh
    [[nodiscard]] size_t numNodes() const noexcept { return nodes_.size(); }

    /// Number of faces in the mesh
    [[nodiscard]] size_t numFaces() const noexcept { return faces_.size(); }

    /// Number of cells in the mesh
    [[nodiscard]] size_t numCells() const noexcept { return cells_.size(); }

    /// Cell count at startup (used by CellData/FaceData)
    [[nodiscard]] static size_t cellCount() noexcept { return cellCount_; }

    /// Face count at startup (used by CellData/FaceData)
    [[nodiscard]] static size_t faceCount() noexcept { return faceCount_; }

private:

    /// Mesh node coordinates
    std::vector<Vector> nodes_;

    /// Mesh face topology
    std::vector<Face> faces_;

    /// Mesh cell topology
    std::vector<Cell> cells_;

    /// Boundary patch descriptors
    std::vector<BoundaryPatch> patches_;

    /// Cell count for field container initialization
    static inline size_t cellCount_ = 0;

    /// Face count for field container initialization
    static inline size_t faceCount_ = 0;
};
