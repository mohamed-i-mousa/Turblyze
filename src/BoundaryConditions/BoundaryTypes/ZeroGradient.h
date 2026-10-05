/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ZeroGradient.h
 * @brief Zero-normal-gradient boundary condition
 *
 * @details The face value is the owner-cell value and no diffusive flux
 * a = 1, b = 0; c = 0, d = 0.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "BoundaryType.h"

// **************************** class ZeroGradient ****************************

class ZeroGradient : public BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Default constructor
    ZeroGradient() noexcept = default;

// ***************************** Override Methods *****************************

    /// The boundary condition type name
    [[nodiscard]] std::string_view typeName() const noexcept override
    {
        return "zeroGradient";
    }

    /// Update patch-local linearized boundary coefficients
    void updateCoeffs
    (
        const Mesh& mesh,
        const BoundaryPatch& patch
    ) override;
};
