/******************************************************************************
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
