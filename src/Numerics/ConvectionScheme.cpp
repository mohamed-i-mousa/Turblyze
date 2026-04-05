/******************************************************************************
 * @file ConvectionScheme.cpp
 * @brief Implementation of convection discretization schemes
 *****************************************************************************/

#include "ConvectionScheme.hpp"

#include "LinearInterpolation.hpp"


// ************************* Central Difference *************************

Scalar CentralDifferenceScheme::calculateCorrection
(
    const Face& face,
    const ScalarField& phi,
    const Vector& /*gradPhiP*/,
    const Vector& /*gradPhiN*/,
    Scalar flowRate
) const
{
    Scalar phiFaceCentral = interpolateToFace(face, phi);

    size_t upwindCell =
        (flowRate >= S(0.0)) ? face.ownerCell() : face.neighborCell().value();

    Scalar phiFaceUDS = phi[upwindCell];

    // Deferred correction: flowRate × (φ_central - φ_upwind)
    return flowRate * (phiFaceCentral - phiFaceUDS);
}


// ************************* Second-Order Upwind *************************

Scalar SecondOrderUpwindScheme::calculateCorrection
(
    const Face& face,
    const ScalarField& /*phi*/,
    const Vector& gradPhiP,
    const Vector& gradPhiN,
    Scalar flowRate
) const
{
    // Deferred correction: flowRate × ∇φ_upwind · d_{upwind→face}
    Scalar correction =
        (flowRate >= S(0.0))
      ? dot(gradPhiP, face.dPf())
      : dot(gradPhiN, face.dNf().value());

    return flowRate * correction;
}


// ********************************** Upwind **********************************

Scalar UpwindScheme::calculateCorrection
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
