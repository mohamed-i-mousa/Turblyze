/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PotentialFlow.h
 * @brief Potential flow initialization via Poisson equation
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Initializer.h"
#include "Vector.h"
#include "Scalar.h"

// *************************** Forward Declarations ***************************

class GradientScheme;
class LinearSolver;

// ************************* namespace Initialization *************************

namespace Initialization
{

// **************************** class PotentialFlow ***************************

class PotentialFlow final : public Initializer
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    PotentialFlow
    (
        const Vector& Uinf,
        Scalar pinf,
        Scalar rho,
        const GradientScheme& gradScheme,
        Count numNonOrthoCorr,
        LinearSolver& solver
    ) noexcept
    :
        Uinf_{Uinf},
        pinf_{pinf},
        rho_{rho},
        gradScheme_{gradScheme},
        numNonOrthoCorr_{numNonOrthoCorr},
        solver_{solver}
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

    /// Freestream velocity vector
    Vector Uinf_;

    /// Reference pressure
    Scalar pinf_;

    /// Fluid density
    Scalar rho_;

    /// Gradient reconstruction scheme
    const GradientScheme& gradScheme_;

    /// Non-orthogonal corrector sub-iterations
    Count numNonOrthoCorr_;

    /// Pressure linear solver reference
    LinearSolver& solver_;
};

} // namespace Initialization
