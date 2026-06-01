/******************************************************************************
 * @file SecondOrderUpwindScheme.h
 * @brief Second-order upwind convection scheme
 *
 * @details Second-order upwind adds a deferred correction based on the upwind
 * cell gradient reconstructed from the upwind cell center to the face.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "ConvectionSchemes.h"

// *********************** class SecondOrderUpwindScheme **********************

class SecondOrderUpwindScheme final : public ConvectionSchemes
{
public:

// ****************************** Public Methods ******************************

    [[nodiscard]] Scalar correction
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};
