/******************************************************************************
 * @file CFDApplication.h
 * @brief Top-level application driver for the CFD solver
 *
 * @details This header defines the CFDApplication class, which manages the
 * entire simulation: case file parsing, mesh preparation, boundary condition
 * setup, solver configuration, solution, and output.
 *
 * @class CFDApplication
 *
 * The CFDApplication class provides:
 * - Phase-based simulation workflow (load, mesh, BCs, solve, export)
 * - Ownership of all simulation data (mesh, fields, solver)
 * - Single entry point via run() method
 *****************************************************************************/

#pragma once

#include <string>
#include <string_view>
#include <memory>

#include "Scalar.h"
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "ConvectionScheme.h"

class GradientScheme;
class SIMPLE;
class CaseReader;
class LinearSolver;

class CFDApplication
{
public:

    /// Constructor for CFDApplication
    explicit CFDApplication(const std::string& caseFilePath);

    /// Copy constructor and assignment - Not copyable (unique_ptr members)
    CFDApplication(const CFDApplication&) = delete;
    CFDApplication& operator=(const CFDApplication&) = delete;

    /// Move constructor and assignment
    CFDApplication(CFDApplication&&) noexcept = default;
    CFDApplication& operator=(CFDApplication&&) noexcept = default;

    /// Destructor
    ~CFDApplication() noexcept;

    /// Run the full simulation
    void run();

private:

// Phase methods

    /// Parse all configuration from case file
    void loadCase();
    
    /// Initialize Eigen threading and report OpenMP thread count
    void initParallelism();

    /// Read mesh, compute geometry, correct normals, check quality
    void prepareMesh();

    /// Register boundary conditions from case file
    void setupBoundaryConditions();

    /// Create and configure the SIMPLE solver
    void configureSolver();

    /// Execute the SIMPLE algorithm
    void solve();

    /// Extract fields, compute statistics, print summary
    void postProcess();

    /// Prepare field maps and write VTK output
    void exportResults();

// Helper methods

    /// Create a convection scheme by scheme name
    static std::unique_ptr<ConvectionSchemes>
    createConvectionScheme(const std::string& name);

    /// Create a linear solver
    static std::unique_ptr<LinearSolver> createLinearSolver
    (
        std::string_view name,
        Scalar tolerance,
        int maxIterations
    );

    /// Parse convection schemes from the case file
    ConvectionScheme parseConvectionSchemes() const;

    /// Read per-equation linear-solver settings and validate them
    static void readAndValidateSolverConfig
    (
        const CaseReader& solvers,
        const std::string& key,
        std::string& solverName,
        std::string& preconditioner,
        Scalar& tolerance,
        int& maxIterations
    );

    /// Report "unknown BC type" error
    [[noreturn]] static void unknownBCType
    (
        std::string_view bcType,
        std::string_view fieldName,
        std::string_view patchName,
        std::string_view validList
    );

    /// Verify k/omega/nut wall functions are set as a complete set
    void validateWallFunctionSetup() const;

// Private members

    /// Path to case file
    std::string caseFilePath_;

    /// Case file reader
    std::unique_ptr<CaseReader> caseReader_;

    /// Mesh data
    Mesh mesh_{};

    /// Physical properties
    Scalar rho_ = 0;
    Scalar mu_ = 0;

    /// Initial conditions
    Vector initialVelocity_;
    Scalar initialPressure_ = 0;
    Scalar initialK_ = 0;
    Scalar initialOmega_ = 0;

    /// SIMPLE algorithm parameters
    int maxIterations_ = 0;
    Scalar convergenceTolerance_ = 0;
    Scalar alphaU_ = 0;
    Scalar alphaP_ = 0;
    Scalar alphaK_ = 0;
    Scalar alphaOmega_ = 0;

    /// Turbulence configuration
    bool turbulenceEnabled_ = false;
    std::string turbulenceModel_;
    Scalar turbIntensity_ = S(0.05);
    Scalar hydrDiameter_ = S(0.01);

    /// Mesh configuration
    bool checkQuality_ = true;

    /// OpenMP thread count (must be > 0; defaults to 1 = serial)
    int numThreads_ = 1;

    /// Output configuration
    std::string vtkOutputFilename_;

    /// Enable verbose console output
    bool debug_ = false;

    /// Velocity constraint configuration
    bool velocityConstraintEnabled_ = false;
    Scalar maxVelocityConstraint_ = 0;

    /// Pressure constraint configuration
    bool pressureConstraintEnabled_ = false;
    Scalar minPressureConstraint_ = 0;
    Scalar maxPressureConstraint_ = 0;

    /// Solver components
    BoundaryConditions bcManager_;
    std::unique_ptr<GradientScheme> gradScheme_;
    ConvectionScheme convectionSchemes_;
    std::unique_ptr<SIMPLE> solver_;
};
