/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FixedGradient.h
 * @brief Fixed-gradient (Neumann) boundary condition
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "BoundaryType.h"

// **************************** class FixedGradient ***************************

class FixedGradient final : public BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Construct with the prescribed boundary-normal gradient
    explicit FixedGradient(Scalar gradient) noexcept
    :
        gradient_{gradient}
    {}

// ***************************** Override Methods *****************************

    /// The boundary condition type name
    [[nodiscard]] std::string_view typeName() const noexcept override
    {
        return "fixedGradient";
    }

    /// The boundary flux receives the explicit p'-gradient correction
    [[nodiscard]] bool correctsBoundaryFlux() const noexcept override
    {
        return true;
    }

    /// Prescribed boundary-normal gradient
    [[nodiscard]] Scalar gradient() const noexcept
    {
        return gradient_;
    }

    /// Update patch-local linearized boundary coefficients
    void updateCoeffs
    (
        const Mesh& mesh,
        const BoundaryPatch& patch
    ) override;

// ****************************** Private Members *****************************

private:

    /// Prescribed boundary-normal gradient
    Scalar gradient_;
};
