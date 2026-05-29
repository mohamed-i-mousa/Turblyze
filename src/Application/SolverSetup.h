/******************************************************************************
 * @file SolverAssembly.h
 * @brief Runtime service ownership and SIMPLE solver assembly
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <memory>

// *************************** Forward Declarations ***************************

class BoundaryConditions;
class ConvectionSchemes;
class GradientScheme;
class kOmegaSST;
class LinearSolver;
class Mesh;
class SIMPLE;
struct CaseConfiguration;

// *************************** struct SolverModules ***************************

struct SolverModules
{
    /// Constructor
    SolverModules();

    /// Copy constructor and assignment - Not copyable (unique_ptr members)
    SolverModules(const SolverModules&) = delete;
    SolverModules& operator=(const SolverModules&) = delete;

    /// Move constructor and assignment - Not movable (borrowed references)
    SolverModules(SolverModules&&) = delete;
    SolverModules& operator=(SolverModules&&) = delete;

    /// Destructor
    ~SolverModules() noexcept;

    /// Gradient scheme owned for SIMPLE's borrowed reference
    std::unique_ptr<GradientScheme> gradScheme;

    /// Default convection scheme fallback
    std::unique_ptr<ConvectionSchemes> defaultConvectionScheme;

    /// Momentum equation convection scheme override
    std::unique_ptr<ConvectionSchemes> momentumConvectionScheme;

    /// k equation convection scheme override
    std::unique_ptr<ConvectionSchemes> kConvectionScheme;

    /// omega equation convection scheme override
    std::unique_ptr<ConvectionSchemes> omegaConvectionScheme;

    /// Momentum linear solver
    std::unique_ptr<LinearSolver> momentumSolver;

    /// Pressure-correction linear solver
    std::unique_ptr<LinearSolver> pressureSolver;

    /// k equation linear solver
    std::unique_ptr<LinearSolver> kSolver;

    /// omega equation linear solver
    std::unique_ptr<LinearSolver> omegaSolver;

    /// Turbulence model (nullptr for laminar runs); destroyed after SIMPLE
    std::unique_ptr<kOmegaSST> turbulenceModel;

    /// SIMPLE solver; declared last so it is destroyed first
    std::unique_ptr<SIMPLE> solver;
};

// *************************** namespace SolverSetup **************************

namespace SolverSetup
{

/// Configure runtime services and construct SIMPLE
void configure
(
    SolverModules& modules,
    const Mesh& mesh,
    const BoundaryConditions& boundaryConditions,
    const CaseConfiguration& config
);

/// Print the SIMPLE solver setup banner
void logSetup
(
    const SolverModules& modules,
    const CaseConfiguration& config
);

} // namespace SolverSetup
