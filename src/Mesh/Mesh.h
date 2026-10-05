/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Mesh.h
 * @brief Owning mesh data
 *
 * @details This header defines the Mesh class, which owns all mesh data
 * (nodes, faces, cells, and boundary patches) and provides list-ref views
 * to consumers. It also tracks the partition halo cell count (0 for serial)
 * to distinguish domain cells from halo cells.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <utility>

// Project headers
#include "Integer.h"
#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"
#include "MeshContainers.h"
#include "ErrorHandler.h"

// ******************************** class Mesh ********************************

class Mesh
{
public:

// ************************* Special Member Functions *************************

    /// Default constructor (empty mesh)
    Mesh() = default;

    /// Constructor
    Mesh
    (
        NodeList nodes,
        FaceList faces,
        CellList cells,
        PatchList patches,
        Count numHaloCells = 0
    )
    :
        nodes_(std::move(nodes)),
        faces_(std::move(faces)),
        cells_(std::move(cells)),
        patches_(std::move(patches)),
        numHaloCells_(numHaloCells)
    {}

    /// Copy constructor and assignment - Not copyable (copy is expensive)
    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    /// Move constructor and assignment
    Mesh(Mesh&&) noexcept = default;
    Mesh& operator=(Mesh&&) noexcept = default;

    /// Destructor
    ~Mesh() noexcept = default;

// ***************************** Accessor Methods *****************************

    /// Node coordinates
    [[nodiscard]] const NodeList& nodes() const noexcept
    {
        return nodes_;
    }

    /// Faces in the mesh
    [[nodiscard]] const FaceList& faces() const noexcept
    {
        return faces_;
    }

    /// Cells in the mesh
    [[nodiscard]] const CellList& cells() const noexcept
    {
        return cells_;
    }

    /// Boundary patches in the mesh
    [[nodiscard]] const PatchList& patches() const noexcept
    {
        return patches_;
    }

// *************************** Size Accessor Methods **************************

    /// Number of nodes in the mesh
    [[nodiscard]] Count numNodes() const noexcept
    {
        return nodes_.size();
    }

    /// Number of faces in the mesh
    [[nodiscard]] Count numFaces() const noexcept
    {
        return faces_.size();
    }

    /// Number of cells in the mesh (halo cells included)
    [[nodiscard]] Count numCells() const noexcept
    {
        return cells_.size();
    }

    /// Number of halo cells
    [[nodiscard]] Count numHaloCells() const noexcept
    {
        return numHaloCells_;
    }

    /// Number of domain (solved) cells
    [[nodiscard]] Count numDomainCells() const noexcept
    {
        return cells_.size() - numHaloCells_;
    }

    /// Whether this mesh is decomposed with halo cells
    [[nodiscard]] bool isDecomposed() const noexcept
    {
        return numHaloCells_ > 0;
    }

// ************************** Geometric Query Methods *************************

    /// Distance vector from owner cell center to face center
    [[nodiscard]] Vector dPf(const Face& f) const noexcept
    {
        return f.centroid() - cells_[f.ownerCell()].centroid();
    }

    /// Distance vector from neighbor cell center to face center
    [[nodiscard]] Vector dNf(const Face& f) const
    {
        return 
            f.centroid()
          - cells_[f.neighborCell().value()].centroid();
    }

    /// Distance vector from owner cell center to neighbor cell center
    [[nodiscard]] Vector dPN(const Face& f) const
    {
        return 
            cells_[f.neighborCell().value()].centroid()
          - cells_[f.ownerCell()].centroid();
    }

    /// gDiff = |Ef| / (|Sf| * |dPf|)
    [[nodiscard]] Scalar gDiff(const Face& f) const noexcept
    {
        const Vector Sf = f.normal() * f.projectedArea();
        const Vector dPfVec = dPf(f);
        const Scalar dPfMag = magnitude(dPfVec);
        const Vector ePf = dPfVec / (dPfMag + vSmallValue);
        const Vector Ef = (dot(Sf, Sf) / dot(Sf, ePf)) * ePf;
        return magnitude(Ef) / (f.projectedArea() * (dPfMag + vSmallValue));
    }

// ****************************** Private Members *****************************

private:

    /// Mesh node coordinates
    NodeList nodes_;

    /// Mesh faces
    FaceList faces_;

    /// Mesh cells
    CellList cells_;

    /// Boundary patches
    PatchList patches_;

    /// Number of halo cells in Parallel runs
    Count numHaloCells_ = 0;
};
