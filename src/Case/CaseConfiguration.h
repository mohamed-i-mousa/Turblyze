/******************************************************************************
 * @file CaseConfiguration.h
 * @brief Typed runtime configuration parsed from a case file
 *
 * @details CaseConfiguration stores non-boundary-condition case input in
 * simple POD-style structs. Members are left bare: CaseConfig::loadConfiguration
 * is the single source of truth for every default and assigns each field on
 * every path. Boundary conditions remain patch/field maps and are streamed
 * directly from CaseReader into BoundaryConditions.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Project headers
#include "Integer.h"
#include "Scalar.h"
#include "Vector.h"
#include "StringTypes.h"

// *************************** Forward Declarations ***************************

class CaseReader;

// **************************** struct SchemeConfig ***************************

struct SchemeConfig
{
    /// Gradient reconstruction scheme name (e.g. "leastSquares")
    Name gradientScheme;

    /// Default convection scheme name
    Name defaultScheme;

    /// Momentum convection scheme name; empty means use default
    Name momentumScheme;

    /// k equation convection scheme name; empty means use default
    Name kScheme;

    /// omega equation convection scheme name; empty means use default
    Name omegaScheme;
};

// ************************ struct LinearSolverSettings ***********************

struct LinearSolverSettings
{
    /// Linear solver runtime-selection name
    Name solver;

    /// Preconditioner name from the case file
    Name preconditioner;

    /// Relative residual tolerance
    Scalar tolerance;

    /// Maximum linear-solver iterations
    Count maxIter;
};

// ************************* struct LinearSolverConfig ************************

struct LinearSolverConfig
{
    /// Momentum equation solver settings
    LinearSolverSettings momentum;

    /// Pressure-correction equation solver settings
    LinearSolverSettings pressure;

    /// Turbulent kinetic energy equation solver settings
    LinearSolverSettings k;

    /// Specific dissipation rate equation solver settings
    LinearSolverSettings omega;
};

// ************************* struct CaseConfiguration *************************

struct CaseConfiguration
{
    /// Mesh file path
    FilePath meshFile;

    /// Whether mesh quality checks should run
    bool checkQuality;

    /// OpenMP/Eigen thread count
    Count numThreads;

    /// Fluid density
    Scalar rho;

    /// Dynamic viscosity
    Scalar mu;

    /// Initial velocity
    Vector initialVelocity;

    /// Initial pressure
    Scalar initialPressure;

    /// Initial turbulent kinetic energy
    Scalar initialK;

    /// Initial specific dissipation rate
    Scalar initialOmega;

    /// Maximum SIMPLE iterations
    Count maxIterations;

    /// SIMPLE convergence tolerance
    Scalar convergenceTolerance;

    /// Velocity relaxation factor
    Scalar alphaU;

    /// Pressure relaxation factor
    Scalar alphaP;

    /// k relaxation factor
    Scalar alphaK;

    /// omega relaxation factor
    Scalar alphaOmega;

    /// Whether turbulence is enabled
    bool turbulenceEnabled;

    /// Turbulence model name
    Name turbulenceModel;

    /// Turbulence intensity for calculated inlet/default values
    Scalar turbulenceIntensity;

    /// Hydraulic diameter for calculated inlet/default values
    Scalar hydraulicDiameter;

    /// Enable velocity field constraints
    bool velocityConstraintEnabled;

    /// Enable pressure field constraints
    bool pressureConstraintEnabled;

    /// Maximum allowed velocity magnitude
    Scalar maxVelocityMagnitude;

    /// Minimum allowed pressure
    Scalar minPressure;

    /// Maximum allowed pressure
    Scalar maxPressure;

    /// Enable aerodynamic force calculation on a wall patch
    bool forcesEnabled;

    /// Wall patch to integrate aerodynamic loads over
    Name forcesPatch;

    /// Unit drag direction (typically the freestream direction)
    Vector dragDirection;

    /// Unit lift direction
    Vector liftDirection;

    /// Reference velocity vector for force coefficients
    Vector referenceVelocity;

    /// Reference (frontal) area for force coefficients
    Scalar referenceArea;

    /// VTK output filename
    FilePath vtkOutputFilename;

    /// Enable verbose output
    bool debug;

    /// Convection scheme names
    SchemeConfig schemes;

    /// Linear solver settings
    LinearSolverConfig linearSolvers;
};

// *************************** namespace CaseConfig ***************************

namespace CaseConfig
{

/// Parse and validate all non-boundary-condition case configuration
[[nodiscard]] CaseConfiguration loadConfiguration(const CaseReader& reader);

} // namespace CaseConfig
