/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file SIMPLE.cpp
 * @brief SIMPLE outer iteration: predictor and corrector
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "SIMPLE.h"

// ************************* Special Member Functions *************************

SIMPLE::SIMPLE
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
)
:
    Segregated
    {
        mesh,
        bc,
        timeScheme,
        gradScheme,
        momentumConvectionScheme,
        momentumSolver,
        pressureSolver,
        turbulence,
        initializer,
        deltaT,
        rho,
        mu,
        alphaU,
        alphaP,
        maxIterations,
        convergenceTolerance,
        nNonOrthogonalCorrectors,
        nOuterCorrectors,
        debug
    }
{}

// ******************************* SIMPLE Solve *******************************

bool SIMPLE::outerIteration
(
    const TransientFields* prevStep
)
{
    // Momentum predictor
    solveMomentum(prevStep);

    // Face mass flux with Rhie-Chow interpolation
    updateRhieChowFlowRate(prevStep);

    // Pressure correction p'
    solvePressureCorrection();

    // Correct fields
    correctVelocity();
    correctFlowRate();
    correctPressure();

    // Solve turbulence
    solveTurbulence();

    // Check convergence
    return checkConvergence();
}