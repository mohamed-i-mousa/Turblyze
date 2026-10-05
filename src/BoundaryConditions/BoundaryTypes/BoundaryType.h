/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryType.h
 * @brief Abstract base class for per-(patch, field) boundary conditions
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <string_view>

// Project headers
#include "BoundaryCoeffs.h"
#include "CellData.h"
#include "Scalar.h"

// *************************** Forward Declarations ***************************

class BoundaryPatch;
class Mesh;

// **************************** class BoundaryType ****************************

class BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Copy constructor and assignment - Not copyable (polymorphic identity)
    BoundaryType(const BoundaryType&) = delete;
    BoundaryType& operator=(const BoundaryType&) = delete;

    /// Move constructor and assignment - Not movable (polymorphic identity)
    BoundaryType(BoundaryType&&) = delete;
    BoundaryType& operator=(BoundaryType&&) = delete;

    /// Destructor
    virtual ~BoundaryType() noexcept = default;

// ****************************** Public Methods ******************************

    /// The boundary condition type name
    [[nodiscard]] virtual std::string_view typeName() const noexcept = 0;

    /// Update patch linearized boundary coefficients
    virtual void updateCoeffs
    (
        const Mesh& mesh,
        const BoundaryPatch& patch
    ) = 0;

    /// Re-evaluate the coefficients that depend on the current velocity
    virtual void refreshCoeffs
    (
        const Mesh& /* mesh */,
        const BoundaryPatch& /* patch */,
        const ScalarField& /* Ux */,
        const ScalarField& /* Uy */,
        const ScalarField& /* Uz */
    )
    {}

    /// Linearization coefficients of a patch face
    [[nodiscard]] Scalar a(Index localIdx) const noexcept
    {
        return coeffs_.a[localIdx];
    }

    [[nodiscard]] Scalar b(Index localIdx) const noexcept
    {
        return coeffs_.b[localIdx];
    }

    [[nodiscard]] Scalar c(Index localIdx) const noexcept
    {
        return coeffs_.c[localIdx];
    }

    [[nodiscard]] Scalar d(Index localIdx) const noexcept
    {
        return coeffs_.d[localIdx];
    }

    /// Reconstruct boundary face value: phi_f = a * adjacentValue + b
    [[nodiscard]] Scalar faceValue
    (
        Index localIdx,
        Scalar adjacentValue
    ) const noexcept
    {
        return coeffs_.a[localIdx] * adjacentValue + coeffs_.b[localIdx];
    }

// ****************************** Capability Flags ****************************

    /// Dirichlet-like: the face value is prescribed
    [[nodiscard]] virtual bool fixesValue() const noexcept
    {
        return false;
    }

    /// Mirror plane: no mass flux, and the velocity components couple
    [[nodiscard]] virtual bool isSymmetry() const noexcept
    {
        return false;
    }

    /// The boundary flux receives the explicit p'-gradient correction
    [[nodiscard]] virtual bool correctsBoundaryFlux() const noexcept
    {
        return false;
    }

    /// Wall-treatment marker: the turbulence model owns the physics
    [[nodiscard]] virtual bool isWallModelled() const noexcept
    {
        return false;
    }

    /// Whether this boundary condition represents a physical wall
    [[nodiscard]] virtual bool isWall() const noexcept
    {
        return false;
    }

// ***************************** Protected Methods ****************************

protected:

    /// Default constructor for derived classes
    BoundaryType() noexcept = default;

// ***************************** Protected Members ****************************

    /// Patch linearized boundary coefficients
    BoundaryCoeffs coeffs_;
};
