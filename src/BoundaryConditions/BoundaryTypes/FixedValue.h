/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FixedValue.h
 * @brief Fixed-value (Dirichlet) boundary condition
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "BoundaryType.h"

// ***************************** class FixedValue *****************************

class FixedValue : public BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Construct with the prescribed boundary value
    explicit FixedValue(Scalar value) noexcept
    :
        value_{value}
    {}

// ***************************** Override Methods *****************************

    /// The boundary condition type name
    [[nodiscard]] std::string_view typeName() const noexcept override
    {
        return "fixedValue";
    }

    /// Dirichlet-like: value is fixed
    [[nodiscard]] bool fixesValue() const noexcept override
    {
        return true;
    }

    /// Prescribed boundary value
    [[nodiscard]] Scalar value() const noexcept
    {
        return value_;
    }

    /// Update patch-local linearized boundary coefficients
    void updateCoeffs
    (
        const Mesh& mesh,
        const BoundaryPatch& patch
    ) override;

// ****************************** Private Members *****************************

private:

    /// Prescribed boundary value
    Scalar value_;
};
