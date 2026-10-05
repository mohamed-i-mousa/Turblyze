/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file NoSlip.h
 * @brief No-slip boundary condition, an alias of FixedValue with value 0
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "FixedValue.h"

// ******************************* class NoSlip *******************************

class NoSlip final : public FixedValue
{
public:

// ************************* Special Member Functions *************************

    /// Construct a zero-value velocity boundary condition
    NoSlip() noexcept
    :
        FixedValue{S(0.0)}
    {}

// ***************************** Override Methods *****************************

    /// The boundary condition type name
    [[nodiscard]] std::string_view typeName() const noexcept override
    {
        return "noSlip";
    }

    /// Whether this boundary condition represents a physical wall
    [[nodiscard]] bool isWall() const noexcept override
    {
        return true;
    }
};