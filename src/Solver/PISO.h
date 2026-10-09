/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PISO.h
 * @brief PISO algorithm for transient incompressible Navier-Stokes equations
 *
 * @details PISO (Pressure-Implicit with Splitting of Operators) is a
 * transient segregated algorithm. Each outer iteration runs one SIMPLE
 * predictor step followed by a fixed number of explicit PRIME corrector
 * steps: the momentum equation is re-assembled with the current flux and
 * advanced explicitly (a single Jacobi sweep) before each pressure correction
 * PISO is transient-only and is run without under-relaxation (alpha = 1)
 *
 * @class PISO
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Segregated.h"

// *************************** Forward Declarations ***************************

class Initializer;

// ******************************** class PISO ********************************

class PISO final : public Segregated
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    PISO
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
        Count nPrimeCorrectors,
        bool debug
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    PISO(const PISO&) = delete;
    PISO& operator=(const PISO&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    PISO(PISO&&) = delete;
    PISO& operator=(PISO&&) = delete;

    /// Destructor
    ~PISO() noexcept override = default;

// ****************************** Private Members *****************************

private:

    /// Number of explicit PRIME correctors per outer iteration
    Count nPrimeCorrectors_;

    /// Scratch read buffers for the explicit Jacobi velocity sweep
    ScalarField UxStar_{mesh_};
    ScalarField UyStar_{mesh_};
    ScalarField UzStar_{mesh_};

// ****************************** Private Methods *****************************

    /// One PISO outer iteration
    [[nodiscard]] bool outerIteration
    (
        const TransientFields* prevStep
    ) override;

    /// Algorithm label for banners and convergence messages
    [[nodiscard]] Name algorithmName() const noexcept override
    {
        return "PISO";
    }

    /// Explicit PRIME momentum
    void solveMomentumExplicit(const TransientFields* prevStep);
};