/******************************************************************************
 * @file ConvectionScheme.cpp
 * @brief Implementation of convection discretization schemes
 *****************************************************************************/

#include "ConvectionScheme.hpp"
#include "LinearInterpolation.hpp"


// ****************************** Base Class ******************************

Scalar ConvectionScheme::calculateCorrection
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
    if (face.isBoundary()) return S(0.0);

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
    const ScalarField& phi,
    const Vector& gradPhiP,
    const Vector& gradPhiN,
    Scalar flowRate
) const
{
    if (face.isBoundary()) return S(0.0);

    // Determine upwind cell based on flow direction
    size_t upwindCell =
        (flowRate >= S(0.0))
      ? face.ownerCell()
      : face.neighborCell().value();

    // Gradient-reconstructed face value: φ_upwind + ∇φ · d
    Scalar phiFaceSOU =
        (upwindCell == face.ownerCell())
      ? phi[upwindCell] + dot(gradPhiP, face.dPf())
      : phi[upwindCell] + dot(gradPhiN, face.dNf().value());

    // Deferred correction: flowRate × (φ_SOU - φ_UDS)
    return flowRate * (phiFaceSOU - phi[upwindCell]);
}
