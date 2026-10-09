/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Uniform.h
 * @brief Uniform flow field initialization
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Initializer.h"
#include "Vector.h"
#include "Scalar.h"

// ************************* namespace Initialization *************************

namespace Initialization
{

// ****************************** class Uniform *******************************

class Uniform final : public Initializer
{
public:

// ************************* Special Member Functions *************************

    /// Constructor with uniform velocity and pressure
    Uniform(const Vector& U, Scalar p) noexcept
    :
        U_{U},
        p_{p}
    {}

// *************************** Initialization Method **************************

    void initialize
    (
        const Mesh& mesh,
        BoundaryConditions& bc,
        ScalarField& Ux,
        ScalarField& Uy,
        ScalarField& Uz,
        ScalarField& p
    ) const override;

// ****************************** Private Members *****************************

private:

    /// Uniform initial velocity
    Vector U_;

    /// Uniform initial pressure
    Scalar p_;
};

} // namespace Initialization
