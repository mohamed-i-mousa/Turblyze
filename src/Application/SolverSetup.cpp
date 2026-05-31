/******************************************************************************
 * @file SolverSetup.cpp
 * @brief Runtime service ownership and SIMPLE solver assembly
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "SolverSetup.h"

// Standard library headers
#include <iomanip>
#include <sstream>
#include <string_view>

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
#include "TurbulenceModel.h"
#include "Laminar.h"
#include "kOmegaSST.h"

// ***************************** Internal Helpers *****************************

namespace
{

std::unique_ptr<ConvectionSchemes> makeConvectionScheme
(
    std::string_view name
)
{
    if (name == "Upwind")
    {
        return std::make_unique<UpwindScheme>();
    }

    if (name == "CentralDifference")
    {
        return std::make_unique<CentralDifferenceScheme>();
    }

    if (name == "SecondOrderUpwind")
    {
        return std::make_unique<SecondOrderUpwindScheme>();
    }

    FatalError("Unknown convection scheme: " + std::string(name));
}


void makeConvectionSchemes
(
    SolverModules& modules,
    const SchemeConfig& config
)
{
    modules.defaultConvectionScheme =
        makeConvectionScheme(config.defaultName);
    modules.momentumConvectionScheme.reset();
    modules.kConvectionScheme.reset();
    modules.omegaConvectionScheme.reset();

    if (!config.momentumName.empty())
    {
        modules.momentumConvectionScheme =
            makeConvectionScheme(config.momentumName);
    }
    if (!config.kName.empty())
    {
        modules.kConvectionScheme =
            makeConvectionScheme(config.kName);
    }
    if (!config.omegaName.empty())
    {
        modules.omegaConvectionScheme =
            makeConvectionScheme(config.omegaName);
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
    if (std::string_view{config.solver} == BiCGSTAB::typeName)
    {
        return std::make_unique<BiCGSTAB>
        (
            config.tolerance,
            config.maxIter
        );
    }

    if (std::string_view{config.solver} == PCG::typeName)
    {
        return std::make_unique<PCG>
        (
            config.tolerance,
            config.maxIter
        );
    }

    FatalError
    (
        "Unknown linear solver '" + config.solver
      + "'. Supported solvers: BiCGSTAB, PCG."
    );
}


void logLinearSolver
(
    std::string_view field,
    const LinearSolverSettings& config
)
{
    Logger::linearSolverConfigRow
    (
        field,
        std::string_view{config.solver},
        config.tolerance,
        config.maxIter
    );
}

} // namespace

// ************************* Special Member Functions *************************

SolverModules::SolverModules() = default;

SolverModules::~SolverModules() noexcept = default;

// ******************************** Configure *********************************

void SolverSetup::configure
(
    SolverModules& modules,
    const Mesh& mesh,
    const BoundaryConditions& boundaryConditions,
    const CaseConfiguration& config
)
{
    modules.gradScheme =
        std::make_unique<GradientScheme>(mesh, boundaryConditions);

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

    if (config.turbulenceEnabled)
    {
        modules.kSolver =
            makeLinearSolver
            (
                config.linearSolvers.k
            );
        modules.omegaSolver =
            makeLinearSolver
            (
                config.linearSolvers.omega
            );

        modules.turbulenceModel =
            std::make_unique<kOmegaSST>
            (
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
            config.velocityConstraintEnabled,
            config.pressureConstraintEnabled,
            config.maxVelocityMagnitude,
            config.minPressure,
            config.maxPressure,
            config.debug
        );
}

// ********************************* logSetup *********************************

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
    if (config.turbulenceEnabled)
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
    Logger::keyValue("Velocity relaxation", config.alphaU);
    Logger::keyValue("Pressure relaxation", config.alphaP);
    if (config.turbulenceEnabled)
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

    if (config.turbulenceEnabled)
    {
        Logger::subsection("Turbulence initialization");
        Logger::keyValue("Model", std::string_view{config.turbulenceModel});
        Logger::keyValue
        (
            "Wall distance",
            std::string
            {
                modules.turbulenceModel->wallDistanceConverged()
              ? "meshWave converged"
              : "meshWave hit iteration cap (results may be degraded)"
            }
        );
        Logger::keyValue("Fields initialized", std::string{"k, omega, nut"});
    }

    Logger::iterationFooter();
}
