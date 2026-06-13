/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LeastSquares.h
 * @brief Weighted least-squares gradient reconstruction scheme
 *
 * @details This header defines the LeastSquares concrete gradient scheme,
 * which computes cell-centered gradients via a weighted least-squares fit over
 * neighbor cells and boundary faces, using inverse distance squared weighting.
 * The per-cell inverse of the normal-equations matrix (AᵀA) is pre-computed
 * once at construction and cached in compact upper-triangular form.
 *
 * @class LeastSquares
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <vector>
#include <array>

// Project headers
#include "GradientScheme.h"

// **************************** class LeastSquares ****************************

class LeastSquares final : public GradientScheme
{
public:

    using InverseATACache = std::vector<std::array<Scalar, 6>>;

// ************************* Special Member Functions *************************

    /// Construct least-squares scheme and pre-compute inverse ATA
    LeastSquares
    (
        const Mesh& mesh,
        const BoundaryConditions& bc
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    LeastSquares(const LeastSquares&) = delete;
    LeastSquares& operator=(const LeastSquares&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    LeastSquares(LeastSquares&&) = delete;
    LeastSquares& operator=(LeastSquares&&) = delete;

    /// Destructor
    ~LeastSquares() noexcept override = default;

// ***************************** Override Methods *****************************

    /// Calculate gradient at a single cell using weighted least-squares
    [[nodiscard]] Vector cellGradient
    (
        Field field,
        const ScalarField& phi,
        Index cellIdx
    ) const override;

// ************************ Private Methods and Members ***********************

private:

    /// Pre-compute inverse of ATA matrix for each cell
    void precomputeInverseATA();

    /// Cached inverse of ATA per cell {xx, xy, xz, yy, yz, zz}
    InverseATACache invATA_;
};
