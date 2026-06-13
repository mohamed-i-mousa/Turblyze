/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file UpwindScheme.cpp
 * @brief Implementation of the first-order upwind convection scheme
 *****************************************************************************/

// ********************************** Headers *********************************

#include "UpwindScheme.h"

// ****************************** Public Methods ******************************

Scalar UpwindScheme::correction
(
    const Face& /*face*/,
    const ScalarField& /*phi*/,
    const Vector& /*gradPhiP*/,
    const Vector& /*gradPhiN*/,
    Scalar /*flowRate*/
) const
{
    return S(0.0);
}
