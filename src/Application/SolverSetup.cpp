/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file SolverSetup.cpp
 * @brief Runtime service ownership and SIMPLE solver assembly
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "SolverSetup.h"

// Standard library headers
#include <format>

// Project headers
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "TimeScheme.h"
#include "GradientScheme.h"
#include "ConvectionScheme.h"
#include "LinearSolvers.h"
#include "PETScRuntime.h"
#include "MomentumTransport.h"
#include "TurbulenceModel.h"
#include "Laminar.h"
#include "CaseConfiguration.h"
#include "ErrorHandler.h"
#include "Logger.h"
#include "StringTypes.h"
#include "Initializer.h"

// ***************************** Internal Helpers *****************************

namespace
{

void makeConvectionSchemes
(
    SolverModules& modules,
    const SchemeConfig& config
)
{
    modules.defaultConvectionScheme =
        ConvectionScheme::create(config.defaultScheme);

    if (!config.momentumScheme.empty())
    {
        modules.momentumConvectionScheme =
            ConvectionScheme::create(config.momentumScheme);
    }
    if (!config.kScheme.empty())
    {
        modules.kConvectionScheme =
            ConvectionScheme::create(config.kScheme);
    }
    if (!config.omegaScheme.empty())
    {
        modules.omegaConvectionScheme =
            ConvectionScheme::create(config.omegaScheme);
    }
}


const ConvectionScheme& resolveConvectionScheme
(
    const std::unique_ptr<ConvectionScheme>& specific,
    const ConvectionScheme& fallback
)
{
    if (specific)
    {
        return *specific;
    }

    return fallback;
}


std::unique_ptr<LinearSolver> makeLinearSolver
(
    const LinearSolverSettings& config,
    const Name& optionsPrefix
)
{
    return LinearSolver::create
    (
        config.solver,
        config.preconditioner,
        config.tolerance,
        config.maxIter,
        optionsPrefix
    );
}


void logLinearSolver
(
    const Name& fieldName,
    const LinearSolverSettings& config
)
{
    Logger::linearSolverConfigRow
    (
        fieldName,
        config.solver,
        config.preconditioner,
        config.tolerance,
        config.maxIter
    );
}

} // namespace

// ************************* Special Member Functions *************************

SolverModules::SolverModules() = default;

SolverModules::~SolverModules() noexcept = default;

// *************************** namespace SolverSetup **************************

void SolverSetup::configure
(
    SolverModules& modules,
    const Mesh& mesh,
    BoundaryConditions& boundaryConditions,
    const CaseConfiguration& config
)
{
    modules.gradScheme =
        GradientScheme::create
        (
            config.schemes.gradientScheme,
            mesh,
            boundaryConditions
        );

    makeConvectionSchemes(modules, config.schemes);

    modules.timeScheme = TimeScheme::create
    (
        config.time.timeScheme,
        config.time.CrankNicolsonCoeff
    );

    PETScRuntime::insertOptions(config.linearSolvers.petscOptions);

    modules.momentumSolver =
        makeLinearSolver
        (
            config.linearSolvers.momentum,
            "momentum_"
        );
    modules.pressureSolver =
        makeLinearSolver
        (
            config.linearSolvers.pressure,
            "pressure_"
        );

    if (!TurbulenceModel::isLaminar(config.turbulenceModel))
    {
        modules.kSolver = makeLinearSolver(config.linearSolvers.k, "k_");
        modules.omegaSolver =
            makeLinearSolver(config.linearSolvers.omega, "omega_");

        modules.turbulenceModel =
            TurbulenceModel::create
            (
                config.turbulenceModel,
                mesh,
                boundaryConditions,
                *modules.timeScheme,
                *modules.gradScheme,
                resolveConvectionScheme
                (
                    modules.kConvectionScheme,
                    *modules.defaultConvectionScheme
                ),
                *modules.kSolver,
                resolveConvectionScheme
                (
                    modules.omegaConvectionScheme,
                    *modules.defaultConvectionScheme
                ),
                *modules.omegaSolver,
                config.time.timeStep,
                config.mu / config.rho,
                config.initialK,
                config.initialOmega,
                config.alphaK,
                config.alphaOmega,
                config.roughWall,
                config.debug
            );
    }
    else
    {
        modules.turbulenceModel =
            std::make_unique<Laminar>(mesh, boundaryConditions, config.mu / config.rho);
    }

    const auto initializer =
        Initializer::create(config, *modules.gradScheme, *modules.pressureSolver);

    modules.solver =
        MomentumTransport::create
        (
            config.algorithm,
            mesh,
            boundaryConditions,
            *modules.timeScheme,
            *modules.gradScheme,
            resolveConvectionScheme
            (
                modules.momentumConvectionScheme,
                *modules.defaultConvectionScheme
            ),
            *modules.momentumSolver,
            *modules.pressureSolver,
            *modules.turbulenceModel,
            *initializer,
            config.time.timeStep,
            config.rho,
            config.mu,
            config.alphaU,
            config.alphaP,
            config.maxIterations,
            config.convergenceTolerance,
            config.nNonOrthogonalCorrectors,
            config.time.nOuterCorrectors,
            config.nPrimeCorrectors,
            config.debug
        );
}


void SolverSetup::logSetup
(
    const SolverModules& modules,
    const CaseConfiguration& config
)
{
    Logger::sectionHeader(std::format("Initializing {} Solver", config.algorithm));

    Logger::subsection("Flow initialization");
    Logger::keyValue("Mode", Name{config.initializationType});

    Logger::subsection("Linear solvers");
    Logger::linearSolverConfigHeader();
    logLinearSolver
    (
        "U",
        config.linearSolvers.momentum
    );
    logLinearSolver
    (
        "p",
        config.linearSolvers.pressure
    );
    if (!TurbulenceModel::isLaminar(config.turbulenceModel))
    {
        logLinearSolver("k", config.linearSolvers.k);
        logLinearSolver
        (
            "omega",
            config.linearSolvers.omega
        );
    }

    Logger::subsection(std::format("{} controls", config.algorithm));

    if (!modules.timeScheme->isTransient())
    {
        Logger::keyValue("Max iterations", config.maxIterations);
        Logger::keyValue("Convergence tolerance", config.convergenceTolerance);
    }
    Logger::keyValue
    (
        "Non-orth correctors",
        config.nNonOrthogonalCorrectors
    );
    Logger::keyValue("Velocity relaxation", config.alphaU);
    Logger::keyValue("Pressure relaxation", config.alphaP);
    if (!TurbulenceModel::isLaminar(config.turbulenceModel))
    {
        Logger::keyValue("k relaxation", config.alphaK);
        Logger::keyValue("omega relaxation", config.alphaOmega);
    }

    if (modules.timeScheme->isTransient())
    {
        Logger::subsection("Time integration");
        Logger::keyValue("Scheme", Name{config.time.timeScheme});
        Logger::keyValue("Time step", config.time.timeStep, "s");
        Logger::keyValue("Total time", config.time.totalTime, "s");
        Logger::keyValue("Write interval", config.time.writingIntervals);
        Logger::keyValue("Outer correctors", config.time.nOuterCorrectors);
        if (config.algorithm == "PISO")
        {
            Logger::keyValue("PRIME correctors", config.nPrimeCorrectors);
        }
    }

    if (!TurbulenceModel::isLaminar(config.turbulenceModel))
    {
        Logger::subsection("Turbulence initialization");
        Logger::keyValue("Model", Name{config.turbulenceModel});
        Logger::keyValue
        (
            "Wall distance",
            Message
            {
                modules.turbulenceModel->wallDistanceConverged()
              ? "meshWave converged"
              : "meshWave hit iteration cap (results may be degraded)"
            }
        );
        Logger::keyValue("Fields initialized", Message{"k, omega, nut"});
    }

    Logger::iterationFooter();
}
