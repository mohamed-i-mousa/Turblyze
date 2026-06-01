/******************************************************************************
 * @file UpwindScheme.h
 * @brief First-order upwind convection scheme
 *
 * @details The upwind scheme uses the implicit first-order upwind matrix
 * coefficients without adding a deferred correction term.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "ConvectionSchemes.h"

// **************************** class UpwindScheme ****************************

class UpwindScheme final : public ConvectionSchemes
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
