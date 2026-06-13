/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file SIMPLE.h
 * @brief SIMPLE algorithm for incompressible Navier-Stokes equations
 *
 * @details This header implements the Semi-Implicit Method for Pressure-Linked
 * Equations (SIMPLE) for solving incompressible flow on unstructured finite
 * volume meshes. The algorithm handles velocity-pressure coupling through
 * pressure correction and turbulence through an abstract TurbulenceModel
 *
 * @class SIMPLE
 * - Pressure-velocity coupling via SIMPLE algorithm
 * - Rhie-Chow interpolation for collocated grid arrangement
 * - Momentum equations with implicit under-relaxation
 * - Pressure correction with mass conservation enforcement
 * - Turbulence model through a TurbulenceModel reference
 * - Field constraints for numerical stability
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <vector>

// Project headers
#include "Mesh.h"
#include "Vector.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "FaceData.h"
#include "GradientScheme.h"
#include "ConvectionSchemes.h"
#include "Matrix.h"
#include "LinearSolvers.h"
#include "Constraint.h"
#include "StringTypes.h"

// *************************** Forward Declarations ***************************

class TurbulenceModel;

// ******************************* class SIMPLE *******************************

class SIMPLE
{
public:

    using NameRefList = std::vector<NameRef>;

// ************************* Special Member Functions *************************

    /// Constructor
    SIMPLE
    (
        const Mesh& mesh,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionSchemes& momentumConvectionScheme,
        LinearSolver& momentumSolver,
        LinearSolver& pressureSolver,
        TurbulenceModel& turbulence,
        Scalar rho,
        Scalar mu,
        const Vector& initialVelocity,
        Scalar initialPressure,
        Scalar alphaU,
        Scalar alphaP,
        Count maxIterations,
        Scalar convergenceTolerance,
        Count nNonOrthogonalCorrectors,
        bool velocityConstraintEnabled,
        bool pressureConstraintEnabled,
        Scalar maxVelocityMagnitude,
        Scalar minPressure,
        Scalar maxPressure,
        bool debug
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    SIMPLE(const SIMPLE&) = delete;
    SIMPLE& operator=(const SIMPLE&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    SIMPLE(SIMPLE&&) = delete;
    SIMPLE& operator=(SIMPLE&&) = delete;

    /// Destructor
    ~SIMPLE() noexcept = default;

// ******************************* SIMPLE Solve *******************************

    /// Main solve loop for SIMPLE algorithm
    void solve();

// ***************************** Accessor Methods *****************************

    /// Get velocity field
    [[nodiscard]] const ScalarField& Ux() const noexcept { return Ux_; }
    [[nodiscard]] const ScalarField& Uy() const noexcept { return Uy_; }
    [[nodiscard]] const ScalarField& Uz() const noexcept { return Uz_; }

    /// Get pressure field
    [[nodiscard]] const ScalarField& pressure() const noexcept { return p_; }

// ****************************** Private Members *****************************

private:

// Dependencies

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

    /// Reference to BCs
    const BoundaryConditions& bcManager_;

    /// Reference to gradient scheme
    const GradientScheme& gradientScheme_;

    /// Reference to momentum convection scheme
    const ConvectionSchemes& momentumConvectionScheme_;

    /// Linear solver for momentum equations
    LinearSolver& momentumSolver_;

    /// Linear solver for pressure correction equation
    LinearSolver& pressureSolver_;

    /// Turbulence model
    TurbulenceModel& turbulence_;

    /// Matrix constructor and solver object
    Matrix matrixConstruct_;

// Physical properties

    /// Kinematic viscosity
    Scalar nu_;

// Algorithm parameters

    /// Under-relaxation factor for velocity
    Scalar alphaU_;

    /// Under-relaxation factor for pressure
    Scalar alphaP_;

    /// Maximum number of iterations
    Count maxIterations_;

    /// Convergence tolerance
    Scalar tolerance_;

    /// Non-orthogonal corrector sub-iterations for the p' equation
    Count nNonOrthogonalCorrectors_;

    /// Enable verbose console output
    bool debug_;

// Solution fields

    /// Velocity fields
    ScalarField Ux_;
    ScalarField Uy_;
    ScalarField Uz_;

    /// Pressure field
    ScalarField p_;

    /// Field constraint system
    Constraint constraintSystem_;

    /// Pressure correction field
    ScalarField pCorr_;

    /// Velocity from previous iteration
    ScalarField UxPrev_;
    ScalarField UyPrev_;
    ScalarField UzPrev_;

    /// Face velocity field (current iteration)
    FaceData<Scalar> UxAvgf_;
    FaceData<Scalar> UyAvgf_;
    FaceData<Scalar> UzAvgf_;

    /// Face velocity from previous iteration
    FaceData<Scalar> UxAvgPrevf_;
    FaceData<Scalar> UyAvgPrevf_;
    FaceData<Scalar> UzAvgPrevf_;

    /// Mass flux through faces (Rhie-Chow)
    FaceFluxField RhieChowFlowRate_;

    /// Mass flux from previous iteration
    FaceFluxField RhieChowFlowRatePrev_;

    /// Momentum diagonal coefficients
    ScalarField DU_;

    /// Face Momentum diagonal coefficients
    FaceFluxField DUf_;

    /// Flag to allow one-time computation of DU_
    bool DUComputed_ = false;

// Gradient fields

    /// Pressure gradient field
    VectorField gradP_;

    /// Pressure correction gradient field
    VectorField gradPCorr_;

    /// Velocity gradient tensor field
    TensorField gradU_;

    /// Per-component velocity gradients
    VectorField gradUx_;
    VectorField gradUy_;
    VectorField gradUz_;

// Momentum assembly fields

    /// Effective viscosity (laminar + turbulent)
    ScalarField nuEff_;

    /// Effective viscosity at face centres
    FaceData<Scalar> nuEffFace_;

    /// Momentum source terms
    ScalarField UxSource_;
    ScalarField UySource_;
    ScalarField UzSource_;

// Pressure-correction assembly fields

    /// Mass imbalance source for pressure correction
    ScalarField massImbalance_;

// Residual tracking for convergence

    /// Track pressure correction RMS before reset
    Scalar lastPressureCorrectionRMS_ = S(1e9);

    /// Track turbulence field changes between iterations
    NameRefList turbulenceResidualNames_;
    ScalarList lastTurbulenceResiduals_;

    /// First-iteration reference values for scaled residuals
    Scalar massImbalance0_ = S(0.0);
    Scalar velocityResidual0_ = S(0.0);
    Scalar pressureResidual0_ = S(0.0);
    ScalarList turbulenceResidual0_;

// ****************************** Private Methods *****************************

    /// Solve momentum equations for each velocity component
    void solveMomentumEquations();

    /// Update face mass fluxes using Rhie-Chow interpolation
    void updateRhieChowFlowRate();

    /// Assemble and solve the pressure correction Poisson equation
    void solvePressureCorrection();

    /// Apply SIMPLE velocity correction: U = U* - D*gradPCorr
    void correctVelocity();

    /// Update face mass fluxes from pressure correction gradient
    void correctFlowRate();

    /// Update pressure with under-relaxation and reset pCorr
    void correctPressure();

    /// Solve the turbulence transport equations for one iteration
    void solveTurbulence();

    /// Check convergence against scaled residual tolerance
    [[nodiscard]] bool checkConvergence();

    /// Compute limited velocity gradients
    void updateVelocityGradients();

    /// Add Σf (νEff)f · (∇U)f^T · Sf to momentum source terms
    void addTransposeGradientSource();

    /// Solve a single momentum component equation
    void solveMomentumComponent
    (
        TransportEquation& eq,
        const ScalarField& componentPrev
    );

    /// Compute mass imbalance across domain
    [[nodiscard]] Scalar massImbalance() const noexcept;

    /// Compute velocity residual
    [[nodiscard]] Scalar velocityResidual() const noexcept;

    /// Compute pressure residual
    [[nodiscard]] Scalar pressureResidual() const noexcept;
};
