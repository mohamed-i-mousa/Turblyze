/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file SecondOrderUpwindScheme.cpp
 * @brief Implementation of the second-order upwind convection scheme
 *****************************************************************************/

// ********************************** Headers *********************************

#include "SecondOrderUpwindScheme.h"

// ****************************** Public Methods ******************************

Scalar SecondOrderUpwindScheme::correction
(
    const Face& face,
    const ScalarField& /*phi*/,
    const Vector& gradPhiP,
    const Vector& gradPhiN,
    Scalar flowRate
) const
{
    // Deferred correction: flowRate * grad(phi)_upwind dot d_upwind_to_face
    const Scalar correction =
        (flowRate >= S(0.0))
      ? dot(gradPhiP, face.dPf())
      : dot(gradPhiN, face.dNf().value());

    return flowRate * correction;
}
