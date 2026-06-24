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
#include <iomanip>
#include <sstream>

// Project headers
#include "BoundaryConditions.h"
#include "CaseConfiguration.h"
#include "ConvectionSchemes.h"
#include "ErrorHandler.h"
#include "GradientScheme.h"
#include "LinearSolvers.h"
#include "Logger.h"
#include "Mesh.h"
#include "SIMPLE.h"
#include "StringTypes.h"
#include "TurbulenceModel.h"
#include "Laminar.h"

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
        ConvectionSchemes::create(config.defaultScheme);
    modules.momentumConvectionScheme.reset();
    modules.kConvectionScheme.reset();
    modules.omegaConvectionScheme.reset();

    if (!config.momentumScheme.empty())
    {
        modules.momentumConvectionScheme =
            ConvectionSchemes::create(config.momentumScheme);
    }
    if (!config.kScheme.empty())
    {
        modules.kConvectionScheme =
            ConvectionSchemes::create(config.kScheme);
    }
    if (!config.omegaScheme.empty())
    {
        modules.omegaConvectionScheme =
            ConvectionSchemes::create(config.omegaScheme);
    }
}


const ConvectionSchemes& resolveConvectionScheme
(
    const std::unique_ptr<ConvectionSchemes>& specific,
    const std::unique_ptr<ConvectionSchemes>& fallback
)
{
    if (specific)
    {
        return *specific;
    }

    if (!fallback)
    {
        FatalError("Default convection scheme must be set.");
    }

    return *fallback;
}


std::unique_ptr<LinearSolver> makeLinearSolver
(
    const LinearSolverSettings& config
)
{
    return LinearSolver::create
    (
        config.solver,
        config.tolerance,
        config.maxIter
    );
}


void logLinearSolver
(
    NameRef fieldName,
    const LinearSolverSettings& config
)
{
    Logger::linearSolverConfigRow
    (
        fieldName,
        NameRef{config.solver},
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
    const BoundaryConditions& boundaryConditions,
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

    modules.momentumSolver =
        makeLinearSolver
        (
            config.linearSolvers.momentum
        );
    modules.pressureSolver =
        makeLinearSolver
        (
            config.linearSolvers.pressure
        );

    if (!TurbulenceModel::isLaminar(config.turbulenceModel))
    {
        modules.kSolver = makeLinearSolver(config.linearSolvers.k);
        modules.omegaSolver = makeLinearSolver(config.linearSolvers.omega);

        modules.turbulenceModel =
            TurbulenceModel::create
            (
                config.turbulenceModel,
                mesh,
                boundaryConditions,
                *modules.gradScheme,
                resolveConvectionScheme
                (
                    modules.kConvectionScheme,
                    modules.defaultConvectionScheme
                ),
                *modules.kSolver,
                resolveConvectionScheme
                (
                    modules.omegaConvectionScheme,
                    modules.defaultConvectionScheme
                ),
                *modules.omegaSolver,
                config.mu / config.rho,
                config.initialK,
                config.initialOmega,
                config.alphaK,
                config.alphaOmega,
                config.debug
            );
    }
    else
    {
        modules.turbulenceModel =
            std::make_unique<Laminar>(mesh, config.mu / config.rho);
    }

    modules.solver =
        std::make_unique<SIMPLE>
        (
            mesh,
            boundaryConditions,
            *modules.gradScheme,
            resolveConvectionScheme
            (
                modules.momentumConvectionScheme,
                modules.defaultConvectionScheme
            ),
            *modules.momentumSolver,
            *modules.pressureSolver,
            *modules.turbulenceModel,
            config.rho,
            config.mu,
            config.initialVelocity,
            config.initialPressure,
            config.alphaU,
            config.alphaP,
            config.maxIterations,
            config.convergenceTolerance,
            config.nNonOrthogonalCorrectors,
            config.velocityConstraintEnabled,
            config.pressureConstraintEnabled,
            config.maxVelocityMagnitude,
            config.minPressure,
            config.maxPressure,
            config.debug
        );
}


void SolverSetup::logSetup
(
    const SolverModules& modules,
    const CaseConfiguration& config
)
{
    Logger::sectionHeader("Initializing SIMPLE Solver");

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

    Logger::subsection("SIMPLE controls");
    Logger::keyValue("Max iterations", config.maxIterations);
    Logger::keyValue("Convergence tolerance", config.convergenceTolerance);
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

    if
    (
        config.velocityConstraintEnabled
     || config.pressureConstraintEnabled
    )
    {
        Logger::subsection("Field constraints");
        if (config.velocityConstraintEnabled)
        {
            std::ostringstream os;
            os << std::scientific << std::setprecision(6)
               << config.maxVelocityMagnitude << " m/s";
            Logger::keyValue("Velocity max", os.str());
        }
        if (config.pressureConstraintEnabled)
        {
            std::ostringstream os;
            os << "[" << std::scientific << std::setprecision(6)
               << config.minPressure << ", "
               << config.maxPressure << "] Pa";
            Logger::keyValue("Pressure range", os.str());
        }
    }

    if (!TurbulenceModel::isLaminar(config.turbulenceModel))
    {
        Logger::subsection("Turbulence initialization");
        Logger::keyValue("Model", NameRef{config.turbulenceModel});
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
