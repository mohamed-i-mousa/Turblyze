/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CaseConfiguration.h
 * @brief Typed runtime configuration parsed from a case file
 *
 * @details CaseConfiguration stores non-boundary-condition case input in
 * simple structs. Members are left bare:
 * CaseConfig::loadConfiguration is the single source of truth for every
 * default and assigns each field on every path. Boundary conditions remain
 * patch/field maps and are streamed directly from CaseReader into
 * BoundaryConditions.
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

    /// PETSc options-database string, e.g. "-pc_type icc -ksp_view"
    Name petscOptions;
};

// **************************** struct TimeControl ***************************

struct TimeControl
{
    /// Time scheme name
    Name timeScheme;

    /// Time step size [s]
    Scalar timeStep;

    /// Total simulated time [s]
    Scalar totalTime;

    /// Write output every N time steps
    Count writingIntervals;

    /// Fixed number of SIMPLE outer correctors per time step
    Count nOuterCorrectors;

    /// Crank-Nicolson coefficient in [0, 1]
    Scalar CrankNicolsonCoeff;
};

// ************************* struct CaseConfiguration *************************

struct CaseConfiguration
{
    /// Mesh file path
    FilePath meshFile;

    /// Whether mesh quality checks should run
    bool checkQuality;

    /// Velocity-coupling algorithm
    Name algorithm;

    /// Fluid density
    Scalar rho;

    /// Dynamic viscosity
    Scalar mu;

    /// Flow field initialization mode (default: "Uniform")
    Name initializationType = "Uniform";

    /// Potential flow non-orthogonal corrector passes
    Count potentialFlowCorrectors = 3;

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

    /// Non-orthogonal corrector sub-iterations for the p' equation
    Count nNonOrthogonalCorrectors;

    /// Number of explicit PRIME correctors per outer iteration
    Count nPrimeCorrectors;

    /// Velocity relaxation factor
    Scalar alphaU;

    /// Pressure relaxation factor
    Scalar alphaP;

    /// k relaxation factor
    Scalar alphaK;

    /// omega relaxation factor
    Scalar alphaOmega;

    /// Turbulence model name (Laminar or kOmegaSST)
    Name turbulenceModel;

    /// Turbulence intensity for calculated inlet/default values
    Scalar turbulenceIntensity;

    /// Hydraulic diameter for calculated inlet/default values
    Scalar hydraulicDiameter;

    /// Enable SST F3 rough-wall blending
    bool roughWall;

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

    /// Transient control
    TimeControl time;
};

// *************************** namespace CaseConfig ***************************

namespace CaseConfig
{

/// Parse and validate all non-boundary-condition case configuration
[[nodiscard]] CaseConfiguration loadConfiguration(const CaseReader& reader);

} // namespace CaseConfig
