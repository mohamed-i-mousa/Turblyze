/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Matrix.h
 * @brief Matrix assembly for finite volume discretization
 *
 * @details This header defines the Matrix class, which handles assembly of
 * sparse linear systems (Ax=b) for any transport equation type: momentum,
 * pressure correction, and turbulence, described by a TransportEquation
 * bundle (see TransportEquation.h).
 *
 * Coefficients are staged in a COO (coordinate) value array whose sparsity
 * pattern is fixed once at construction from the mesh topology and handed
 * to PETSc via MatSetPreallocationCOO with GLOBAL indices (rank-major
 * owned cells; ghost columns through the mesh's ghost global ids). The
 * slot layout is
 *
 *     [ 2 per internal face | 1 per processor face | 1 per owned cell ]
 *       off-diagonal pairs    local-row off-diags    diagonals
 *
 * A processor face carries only the LOCAL cell's row here — the mirror
 * entry is assembled by the neighbor rank as its own processor face. Both
 * off-diagonal COO writes are still unconditional, so one extra scribble
 * slot at the very end gives the remote side a harmless sink to land in;
 * PETSc never sees it. This rank's matrix rows, RHS, and diagonal cover
 * owned cells only.
 *
 * Assembly is one scatter pass per face group: owned rows are primed with
 * the cell source and any transient term, then every face adds its share
 * to the rows it touches. A face whose far side is a ghost adds nothing
 * there — that row belongs to the neighbor rank.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <span>
#include <vector>

// External library headers
#include <petscmat.h>
#include <petscvec.h>

// Project headers
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "Integer.h"
#include "TransportEquation.h"

// ******************************* class Matrix *******************************

class Matrix
{
public:

// ************************* Special Member Functions *************************

    /// Constructor: builds the PETSc system
    Matrix
    (
        const Mesh& mesh,
        const BoundaryConditions& boundaryConds
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    Matrix(Matrix&&) = delete;
    Matrix& operator=(Matrix&&) = delete;

    /// Destructor releases the PETSc matrix and RHS vector
    ~Matrix() noexcept;

// ****************************** Public Methods ******************************

    /// Build transport equation coefficients
    void buildMatrix(const TransportEquation& equation);

    /// Apply Patankar implicit under-relaxation to the staged values
    void relax(Scalar alpha, const ScalarField& phiPrevIter);

    /// Fix matrix rows to impose known cell values (ghosts included)
    void setValues
    (
        const IndexList& cellIndices,
        const ScalarList& values,
        const ScalarField& ghostFractions,
        const ScalarField& ghostValues,
        const ScalarList& fractions
    );

    /// Push the staged values into the PETSc matrix
    void assemble();

    /// One explicit Jacobi sweep of the staged system
    void explicitJacobiUpdate
    (
        const ScalarField& phiOld,
        ScalarField& phiNew
    ) const;

// ***************************** Accessor Methods *****************************

    /// Get the assembled PETSc matrix
    [[nodiscard]] Mat matrixA() const noexcept
    {
        return matrixA_;
    }

    /// Get the PETSc right-hand side vector
    [[nodiscard]] Vec rhsVec() const noexcept
    {
        return rhsVec_;
    }

    /// Get right-hand side vector b (const)
    [[nodiscard]] const ScalarList& vectorB() const noexcept
    {
        return vectorB_;
    }

    /// Get right-hand side vector b
    [[nodiscard]] ScalarList& vectorB() noexcept
    {
        return vectorB_;
    }

    /// Get the matrix diagonal (const)
    [[nodiscard]] std::span<const Scalar> diagonal() const noexcept
    {
        return {cooValues_.data() + diagOffset_, mesh_.numDomainCells()};
    }

    /// Get the staged matrix diagonal
    [[nodiscard]] std::span<Scalar> diagonal() noexcept
    {
        return {cooValues_.data() + diagOffset_, mesh_.numDomainCells()};
    }

// ****************************** Private Members *****************************

private:

    /// Mesh data reference
    const Mesh& mesh_;

    /// Boundary data references
    const BoundaryConditions& bcManager_;

    /// PETSc sparse matrix
    Mat matrixA_ = nullptr;

    /// PETSc RHS vectorB_ storage
    Vec rhsVec_ = nullptr;

    /// Staged COO values: [internal | processor | diagonals | scribble]
    ScalarList cooValues_;

    /// Right-hand side staging
    ScalarList vectorB_;

    /// Face indices of every internal face, in slot order
    IndexList internalFaces_;

    /// Face indices of every processor (inter-rank cut) face, slot order
    IndexList processorFaces_;

    /// Face indices of every boundary face (no COO slot of their own)
    IndexList boundaryFaces_;

    /// faceIdx -> its COO slot (internal faces: +1 is the reverse)
    IndexList faceSlot_;

    /// First diagonal slot in cooValues_
    Index diagOffset_ = 0;

    /// First scribble slot: absorbs processor faces' remote-side writes
    Index scribbleSlot_ = 0;

    /// Relaxation factor from last relax() call (0 = not relaxed)
    Scalar lastRelaxationFactor_ = S(0.0);

    /// Threshold below which f/(1-f) overflows
    inline static const Scalar rootSmallValue_ = std::sqrt(smallValue);


// ****************************** Private Methods *****************************

private:

    /// Add a two-cell face's contribution to the owned rows it touches
    void assembleInternalFace
    (
        Index cooOwnerSlot,
        Index cooNeighborSlot,
        const Face& face,
        const TransportEquation& equation
    );

    /// Add a boundary face's contribution to its owner cell's row
    void assembleBoundaryFace
    (
        const Face& face,
        const TransportEquation& equation,
        const BoundaryType& bc
    );
};
