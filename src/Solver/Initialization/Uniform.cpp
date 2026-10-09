/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Uniform.cpp
 * @brief Implementation of uniform flow field initialization
 *****************************************************************************/

// ********************************** Headers *********************************

#include "Uniform.h"

// ************************* namespace Initialization *************************

namespace Initialization
{

// *************************** Initialization Method **************************

void Uniform::initialize
(
    const Mesh& /*mesh*/,
    BoundaryConditions& /*bc*/,
    ScalarField& Ux,
    ScalarField& Uy,
    ScalarField& Uz,
    ScalarField& p
) const
{
    Ux.setAll(U_.x());
    Uy.setAll(U_.y());
    Uz.setAll(U_.z());
    p.setAll(p_);
}

} // namespace Initialization
