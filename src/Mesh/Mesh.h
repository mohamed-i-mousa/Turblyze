/******************************************************************************
 * @file Mesh.h
 * @brief Owning mesh data and lightweight view provider
 *
 * @details This header defines the Mesh class, which owns all mesh data
 * (nodes, faces, cells, and boundary patches) and provides list-ref views
 * to consumers. It also manages cell and face counts used by field
 * containers at construction time.
 *
 * @class Mesh
 * - Owns nodes, faces, cells, and boundary patches via mesh list members
 * - Const list-ref accessors for read-only consumers (const Mesh&)
 * - Mutable list-ref accessors for the mesh preparation phase
 * - Static retrieval of cell and face counts
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <utility>

// Project headers
#include "MeshContainers.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"
#include "ErrorHandler.h"
#include "Integer.h"

// ******************************** class Mesh ********************************

class Mesh
{
public:

// ************************* Special Member Functions *************************

    /// Construct empty mesh
    Mesh() = default;

    /// Construct mesh by taking ownership of data vectors
    Mesh
    (
        NodeList nodes,
        FaceList faces,
        CellList cells,
        PatchList patches
    )
    :
        nodes_(std::move(nodes)),
        faces_(std::move(faces)),
        cells_(std::move(cells)),
        patches_(std::move(patches))
    {
        if (cellCount_ != 0 || faceCount_ != 0)
        {
            FatalError
            (
                "Mesh: a populated mesh already exists. "
                "Only one Mesh instance may carry data per program."
            );
        }

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

// ************************** Const Accessor Methods **************************

    /// Read-only view of node coordinates
    [[nodiscard]] NodeListRef nodes() const noexcept
    {
        return nodes_;
    }

    /// Read-only view of faces
    [[nodiscard]] FaceListRef faces() const noexcept
    {
        return faces_;
    }

    /// Read-only view of cells
    [[nodiscard]] CellListRef cells() const noexcept
    {
        return cells_;
    }

    /// Read-only view of boundary patches
    [[nodiscard]] PatchListRef patches() const noexcept
    {
        return patches_;
    }

// ************************* Mutable Accessor Methods *************************

    /// Mutable view of faces
    [[nodiscard]] MutableFaceListRef faces() noexcept
    {
        return faces_;
    }

    /// Mutable view of cells
    [[nodiscard]] MutableCellListRef cells() noexcept
    {
        return cells_;
    }

// *************************** Size Accessor Methods **************************

    /// Number of nodes in the mesh
    [[nodiscard]] Count numNodes() const noexcept { return nodes_.size(); }

    /// Number of faces in the mesh
    [[nodiscard]] Count numFaces() const noexcept { return faces_.size(); }

    /// Number of cells in the mesh
    [[nodiscard]] Count numCells() const noexcept { return cells_.size(); }

    /// Cell count at startup (used by CellData/FaceData)
    [[nodiscard]] static Count cellCount() noexcept { return cellCount_; }

    /// Face count at startup (used by CellData/FaceData)
    [[nodiscard]] static Count faceCount() noexcept { return faceCount_; }

// ****************************** Private Members *****************************

private:

    /// Mesh node coordinates
    NodeList nodes_;

    /// Mesh face topology
    FaceList faces_;

    /// Mesh cell topology
    CellList cells_;

    /// Boundary patch descriptors
    PatchList patches_;

    /// Cell count for field container initialization
    static inline Count cellCount_ = 0;

    /// Face count for field container initialization
    static inline Count faceCount_ = 0;
};
