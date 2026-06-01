/******************************************************************************
 * @file CentralDifferenceScheme.h
 * @brief Central difference convection scheme
 *
 * @details Central difference adds a deferred correction from the implicit
 * upwind face value to the linearly interpolated face value.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "ConvectionSchemes.h"

// *********************** class CentralDifferenceScheme **********************

class CentralDifferenceScheme final : public ConvectionSchemes
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
