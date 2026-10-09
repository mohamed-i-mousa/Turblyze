/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PISO.cpp
 * @brief PISO outer iteration: SIMPLE predictor + explicit PRIME correctors
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "PISO.h"

// Standard library headers
#include <utility>

// Project headers
#include "GradientScheme.h"
#include "HaloExchange.h"

// ************************* Special Member Functions *************************

PISO::PISO
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
    },
    nPrimeCorrectors_{nPrimeCorrectors}
{}

// ******************************** PISO Solve ********************************

bool PISO::outerIteration
(
    const TransientFields* prevStep
)
{
    // SIMPLE predictor
    solveMomentum(prevStep);
    updateRhieChowFlowRate(prevStep);
    solvePressureCorrection();
    correctVelocity();
    correctFlowRate();
    correctPressure();

    // PRIME correctors
    for (Count corrector = 0; corrector < nPrimeCorrectors_; ++corrector)
    {
        // Update Pressure gradient
        gradientScheme().fieldGradient(Field::p, pressure(), gradP());
        Halo::exchange({&gradP()});

        // Explicit momentum predictor
        solveMomentumExplicit(prevStep);
        updateRhieChowFlowRate(prevStep);
        solvePressureCorrection();
        correctVelocity();
        correctFlowRate();
        correctPressure();
    }

    // Solve turbulence
    solveTurbulence();

    // Check convergence
    return checkConvergence();
}


void PISO::solveMomentumExplicit(const TransientFields* prevStep)
{
    assembleMomentum();


    TransportEquation equationUx
    {
        .field          = Field::Ux,
        .phi            = Ux(),
        .transient      =
                    ddtTerm
                    (
                        prevStep,
                        &TransientFields::UxPrevStep,
                        &TransientFields::UxDdtPrevStep
                    ),
        .convection     =
                    ConvectionTerm
                    {
                        faceMassFlux(),
                        momentumConvectionScheme()
                    },
        .GammaFace      = nuEffFace(),
        .source         = UxSource(),
        .gradPhi        = gradUx(),
        .gradScheme     = gradientScheme()
    };
    matrixConstruct().buildMatrix(equationUx);
    diagonalDU(0);
    matrixConstruct().explicitJacobiUpdate(Ux(), UxStar_);

    TransportEquation equationUy
    {
        .field          = Field::Uy,
        .phi            = Uy(),
        .transient      = 
                    ddtTerm
                    (
                        prevStep,
                        &TransientFields::UyPrevStep,
                        &TransientFields::UyDdtPrevStep
                    ),
        .convection     =
                    ConvectionTerm
                    {
                        faceMassFlux(),
                        momentumConvectionScheme()
                    },
        .GammaFace      = nuEffFace(),
        .source         = UySource(),
        .gradPhi        = gradUy(),
        .gradScheme     = gradientScheme()
    };
    matrixConstruct().buildMatrix(equationUy);
    diagonalDU(1);
    matrixConstruct().explicitJacobiUpdate(Uy(), UyStar_);

    TransportEquation equationUz
    {
        .field          = Field::Uz,
        .phi            = Uz(),
        .transient      = 
                    ddtTerm
                    (
                        prevStep,
                        &TransientFields::UzPrevStep,
                        &TransientFields::UzDdtPrevStep
                    ),
        .convection     =
                    ConvectionTerm
                    {
                        faceMassFlux(),
                        momentumConvectionScheme()
                    },
        .GammaFace      = nuEffFace(),
        .source         = UzSource(),
        .gradPhi        = gradUz(),
        .gradScheme     = gradientScheme()
    };
    matrixConstruct().buildMatrix(equationUz);
    diagonalDU(2);
    matrixConstruct().explicitJacobiUpdate(Uz(), UzStar_);

    // Swap the swept values
    std::swap(Ux(), UxStar_);
    std::swap(Uy(), UyStar_);
    std::swap(Uz(), UzStar_);

    // The sweeps wrote owned cells only: refresh U ghosts before reading
    Halo::exchange({&Ux(), &Uy(), &Uz()});

    // Rebuild the face momentum diagonal
    buildFaceDiagonal();
}