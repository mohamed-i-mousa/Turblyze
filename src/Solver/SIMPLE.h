/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file SIMPLE.h
 * @brief SIMPLE algorithm for incompressible Navier-Stokes equations
 *
 * @details SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) is a
 * segregated pressure-correction algorithm. It inherits from MomentumTransport
 * and Segregated. its only contribution is the body of one outer iteration:
 * implicit momentum predictor, Rhie-Chow flux, pressure correction, and the
 * velocity, correction, and the velocity, flux, and pressure corrections
 *
 * @class SIMPLE
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Segregated.h"

// *************************** Forward Declarations ***************************

class Initializer;

// ******************************* class SIMPLE *******************************

class SIMPLE final : public Segregated
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    SIMPLE
    (
        const Mesh& mesh,
        BoundaryConditions& bc,
        const TimeScheme& timeScheme,
        const GradientScheme& gradScheme,
        const ConvectionScheme& momentumConvectionScheme,
        LinearSolver& momentumSolver,
        LinearSolver& pressureSolver,
        TurbulenceModel& turbulence,
        const Initializer& initializer,
        Scalar deltaT,
        Scalar rho,
        Scalar mu,
        Scalar alphaU,
        Scalar alphaP,
        Count maxIterations,
        Scalar convergenceTolerance,
        Count nNonOrthogonalCorrectors,
        Count nOuterCorrectors,
        bool debug
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    SIMPLE(const SIMPLE&) = delete;
    SIMPLE& operator=(const SIMPLE&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    SIMPLE(SIMPLE&&) = delete;
    SIMPLE& operator=(SIMPLE&&) = delete;

    /// Destructor
    ~SIMPLE() noexcept override = default;

// ****************************** Private Methods *****************************

private:

    /// One SIMPLE outer iteration
    [[nodiscard]] bool outerIteration
    (
        const TransientFields* prevStep
    ) override;

    /// Algorithm label for banners and convergence messages
    [[nodiscard]] Name algorithmName() const noexcept override
    {
        return "SIMPLE";
    }
};