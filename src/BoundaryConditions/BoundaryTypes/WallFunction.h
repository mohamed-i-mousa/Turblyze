/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file WallFunction.h
 * @brief Wall-function marker boundary condition for turbulence fields
 *
 * @class WallFunction
 * One class covers the kWallFunction/omegaWallFunction/nutWallFunction
 * flavors, distinguished by the stored parse token. Its assembly numerics
 * are inherited zero-gradient; the wall physics (omega wall values, k
 * production override, nut log-law) lives in the turbulence model, which
 * activates it wherever isWallModelled() is true. A flavor graduates to its
 * own class only when it gains distinct boundary-layer behavior.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Project headers
#include "ZeroGradient.h"

// **************************** class WallFunction ****************************

class WallFunction final : public ZeroGradient
{
public:

// ************************* Special Member Functions *************************

    /// Construct with the flavor's name (e.g. "kWallFunction")
    explicit WallFunction(Name typeName) noexcept
    :
        ZeroGradient{},
        typeName_{std::move(typeName)}
    {}

// ***************************** Override Methods *****************************

    /// The boundary condition type name
    [[nodiscard]] std::string_view typeName() const noexcept override
    {
        return typeName_;
    }

    /// Wall-treatment marker: the turbulence model owns the physics
    [[nodiscard]] bool isWallModelled() const noexcept override
    {
        return true;
    }

// ****************************** Private Members *****************************

private:

    /// The wallFunction flavor name
    Name typeName_;
};
