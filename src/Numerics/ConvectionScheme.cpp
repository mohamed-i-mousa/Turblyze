/******************************************************************************
 * @file ConvectionScheme.cpp
 * @brief Implementation of convection discretization schemes
 *****************************************************************************/

#include <algorithm>

#include "ConvectionScheme.hpp"
#include "LinearInterpolation.hpp"


// ****************************** Base Class ******************************

ConvectionScheme::FluxCoefficients
ConvectionScheme::getFluxCoefficients
(
    Scalar flowRate
)
{
    return
    {
        std::max(flowRate, S(0.0)),
        std::min(flowRate, S(0.0))
    };
}

Scalar ConvectionScheme::calculateCorrection
(
    const Face& /*face*/,
    const ScalarField& /*phi*/,
    const Vector& /*grad_phi_P*/,
    const Vector& /*grad_phi_N*/,
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
    const Vector& /*grad_phi_P*/,
    const Vector& /*grad_phi_N*/,
    Scalar flowRate
) const
{
    if (face.isBoundary()) return S(0.0);

    Scalar phi_face_central = interpolateToFace(face, phi);

    size_t upwind_cell =
        (flowRate >= S(0.0)) ? face.ownerCell() : face.neighborCell().value();

    Scalar phi_face_UDS = phi[upwind_cell];

    // Deferred correction: flowRate × (φ_central - φ_upwind)
    return flowRate * (phi_face_central - phi_face_UDS);
}


// ************************* Second-Order Upwind *************************

Scalar SecondOrderUpwindScheme::calculateCorrection
(
    const Face& face,
    const ScalarField& phi,
    const Vector& grad_phi_P,
    const Vector& grad_phi_N,
    Scalar flowRate
) const
{
    if (face.isBoundary()) return S(0.0);

    // Determine upwind cell based on flow direction
    size_t upwind_cell =
        (flowRate >= S(0.0))
      ? face.ownerCell()
      : face.neighborCell().value();

    // Gradient-reconstructed face value: φ_upwind + ∇φ · d
    Scalar phi_face_SOU =
        (upwind_cell == face.ownerCell())
      ? phi[upwind_cell] + dot(grad_phi_P, face.d_Pf())
      : phi[upwind_cell] + dot(grad_phi_N, face.d_Nf().value());

    // Deferred correction: flowRate × (φ_SOU - φ_UDS)
    return flowRate * (phi_face_SOU - phi[upwind_cell]);
}