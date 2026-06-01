/******************************************************************************
 * @file CentralDifferenceScheme.cpp
 * @brief Implementation of the central difference convection scheme
 *****************************************************************************/

// ********************************** Headers *********************************

#include "CentralDifferenceScheme.h"
#include "LinearInterpolation.h"

// ****************************** Public Methods ******************************

Scalar CentralDifferenceScheme::correction
(
    const Face& face,
    const ScalarField& phi,
    const Vector& /*gradPhiP*/,
    const Vector& /*gradPhiN*/,
    Scalar flowRate
) const
{
    const Scalar phiFaceCentral = interpolateToFace(face, phi);

    const size_t upwindCell =
        (flowRate >= S(0.0)) ? face.ownerCell() : face.neighborCell().value();

    const Scalar phiFaceUDS = phi[upwindCell];

    // Deferred correction: flowRate * (phi_central - phi_upwind)
    return flowRate * (phiFaceCentral - phiFaceUDS);
}
