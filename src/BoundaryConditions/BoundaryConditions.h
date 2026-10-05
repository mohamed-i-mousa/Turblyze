/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryConditions.h
 * @brief Owns boundary conditions for the CFD solver
 *
 * @details Registry of the BoundaryType of every (patch, field) pair, and the
 * single point through which the rest of the solver asks boundary questions:
 * face values, face velocity and mass flux, patch classification, and the
 * velocity-dependent coefficient refresh. Consumers never branch on a
 * concrete boundary condition type.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <map>
#include <memory>

// Project headers
#include "BoundaryType.h"
#include "CellData.h"
#include "Field.h"
#include "Scalar.h"
#include "StringTypes.h"

// *************************** Forward Declarations ***************************

class BoundaryPatch;
class Face;
class Mesh;

// ************************* class BoundaryConditions *************************

class BoundaryConditions
{
public:

    using BCs =
        std::map<Name, std::map<Field, std::unique_ptr<BoundaryType>>>;

// ************************* Special Member Functions *************************

    /// Construct empty boundary conditions
    BoundaryConditions() = default;

    /// Constructor from a pre-built map, evaluating the geometric coefficients
    BoundaryConditions(BCs boundaryConditions, const Mesh& mesh);

    /// Copy constructor and assignment - Not copyable (contains unique_ptr)
    BoundaryConditions(const BoundaryConditions&) = delete;
    BoundaryConditions& operator=(const BoundaryConditions&) = delete;

    /// Move constructor and assignment
    BoundaryConditions(BoundaryConditions&&) noexcept = default;
    BoundaryConditions& operator=(BoundaryConditions&&) noexcept = default;

    /// Destructor
    ~BoundaryConditions() noexcept = default;

// ***************************** Accessor Methods *****************************

    /// Get the boundary condition object for a field on a patch name
    [[nodiscard]] const BoundaryType& boundaryType
    (
        const Name& patchName,
        Field field
    ) const;

    /// Get the boundary condition object for a field on a patch
    [[nodiscard]] const BoundaryType& boundaryType
    (
        const BoundaryPatch& patch,
        Field field
    ) const;

    /// Get the boundary condition object for a field on a boundary face
    [[nodiscard]] const BoundaryType& boundaryType
    (
        const Face& face,
        Field field
    ) const;

    /// Reconstructed boundary face value: phi_f = a * ownerValue + b
    [[nodiscard]] Scalar faceValue
    (
        const Face& face,
        Scalar ownerValue,
        Field field
    ) const noexcept;

// ****************************** Public Methods ******************************

    /// Re-evaluate the velocity-dependent coefficients on every physical patch
    void refresh
    (
        const Mesh& mesh,
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz
    );

    /// Print summary of all boundary conditions
    void printSummary() const;

// ****************************** Private Methods *****************************

private:

    /// Boundary condition of a field on a patch, or nullptr if not set
    [[nodiscard]] const BoundaryType* find
    (
        const Name& patchName,
        Field field
    ) const noexcept;

    /// Abort for a missing boundary condition
    [[noreturn]] static void missing(const Name& patchName, Field field);

// ****************************** Private Members *****************************

    /// Nested map: patch name → field → boundary condition object
    BCs boundaryConditions_;
};
