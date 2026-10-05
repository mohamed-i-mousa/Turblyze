/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Symmetry.h
 * @brief Symmetry-plane boundary condition for a velocity component
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "BoundaryType.h"
#include "Field.h"

// ****************************** class Symmetry ******************************

class Symmetry final : public BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Construct for a velocity component (Ux, Uy or Uz)
    explicit Symmetry(Field component) noexcept
    :
        component_{component}
    {}

// ***************************** Override Methods *****************************

    /// The boundary condition type name
    [[nodiscard]] std::string_view typeName() const noexcept override
    {
        return "symmetry";
    }

    /// Mirror plane: no mass flux, and the velocity components couple
    [[nodiscard]] bool isSymmetry() const noexcept override
    {
        return true;
    }

    /// Update the geometry-only coefficients a and c; b and d are zeroed
    void updateCoeffs
    (
        const Mesh& mesh,
        const BoundaryPatch& patch
    ) override;

    /// Update b and d from the other velocity components
    void refreshCoeffs
    (
        const Mesh& mesh,
        const BoundaryPatch& patch,
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz
    ) override;

// ****************************** Private Members *****************************

private:

    /// The velocity component this symmetry condition applies to
    Field component_;
};
