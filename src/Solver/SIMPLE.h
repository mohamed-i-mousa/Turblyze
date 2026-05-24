/******************************************************************************
 * @file SIMPLE.h
 * @brief SIMPLE algorithm for incompressible Navier-Stokes equations
 *
 * @details This header implements the Semi-Implicit Method for Pressure-Linked
 * Equations (SIMPLE) for solving incompressible flow on unstructured finite
 * volume meshes. The algorithm handles velocity-pressure coupling through
 * pressure correction and includes the k-omega SST turbulence modeling.
 *
 * @class SIMPLE
 *
 * The SIMPLE class provides:
 * - Pressure-velocity coupling via SIMPLE algorithm
 * - Rhie-Chow interpolation for collocated grid arrangement
 * - Momentum equations with implicit under-relaxation
 * - Pressure correction with mass conservation enforcement
 * - k-omega SST turbulence model with wall distance calculation
 * - Field constraints for numerical stability
 *****************************************************************************/

#pragma once

#include <vector>
#include <memory>
#include <utility>

#include "Mesh.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "FaceData.h"
#include "GradientScheme.h"
#include "ConvectionScheme.h"
#include "Matrix.h"
#include "LinearSolvers.h"
#include "ErrorHandler.h"
#include "kOmegaSST.h"
#include "Constraint.h"


class SIMPLE
{
public:

    /// Constructor for SIMPLE solver
    SIMPLE
    (
        const Mesh& mesh,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionScheme& convSchemes
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    SIMPLE(const SIMPLE&) = delete;
    SIMPLE& operator=(const SIMPLE&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    SIMPLE(SIMPLE&&) = delete;
    SIMPLE& operator=(SIMPLE&&) = delete;

    /// Destructor
    ~SIMPLE() noexcept = default;

    /// Main SIMPLE algorithm execution
    void solve();

    /// Initialize solution fields
    void initialize
    (
        const Vector& initialVelocity,
        Scalar initialPressure,
        Scalar initialK = S(1e-6),
        Scalar initialOmega = S(1.0),
        bool enableTurbulence = false
    );

// Accessor methods

    /// Get velocity field
    [[nodiscard]] const ScalarField& Ux() const noexcept { return Ux_; }
    [[nodiscard]] const ScalarField& Uy() const noexcept { return Uy_; }
    [[nodiscard]] const ScalarField& Uz() const noexcept { return Uz_; }

    /// Get pressure field
    [[nodiscard]] const ScalarField& pressure() const noexcept { return p_; }

// Setter methods

    /// Set under-relaxation factors
    void setRelaxationFactors
    (
        Scalar alphaU,
        Scalar alphaP,
        Scalar alphaK = S(0.5),
        Scalar alphaOmega = S(0.5)
    ) noexcept
    {
        alphaU_ = alphaU;
        alphaP_ = alphaP;
        alphaK_ = alphaK;
        alphaOmega_ = alphaOmega;
    }

    /// Set convergence tolerance
    void setConvergenceTolerance(Scalar tol) noexcept { tolerance_ = tol; }

    /// Set maximum number of iterations
    void setMaxIterations(int maxIter) noexcept { maxIterations_ = maxIter; }

    /// Enable or disable verbose console output
    void setDebug(bool d) noexcept { debug_ = d; }

    /// Set physical properties
    void setPhysicalProperties(Scalar rho, Scalar mu) noexcept
    {
        if (rho <= S(0.0))
        {
            FatalError("setPhysicalProperties: density must be positive");
        }

        rho_ = rho;
        mu_  = mu;
        nu_  = mu / rho;
    }

    /// Set linear solver for momentum equations
    void setMomentumSolver(std::unique_ptr<LinearSolver> solver) noexcept
    {
        momentumSolver_ = std::move(solver);
    }

    /// Set linear solver for pressure correction equation
    void setPressureSolver(std::unique_ptr<LinearSolver> solver) noexcept
    {
        pressureSolver_ = std::move(solver);
    }

    /// Set linear solvers for turbulence equations
    void setTurbulenceSolvers
    (
        std::unique_ptr<LinearSolver> kSolver,
        std::unique_ptr<LinearSolver> omegaSolver
    ) noexcept;

// Getter methods

    /// Get constraint system
    [[nodiscard]] Constraint& constraintSystem() noexcept
    {
        return *constraintSystem_;
    }

    /// Get turbulent kinetic energy field
    [[nodiscard]] const ScalarField& turbulentKineticEnergy() const noexcept
    {
        return turbulenceModel_->k();
    }

    /// Get specific dissipation rate field
    [[nodiscard]] const ScalarField& specificDissipationRate() const noexcept
    {
        return turbulenceModel_->omega();
    }

    /// Get turbulent viscosity field
    [[nodiscard]] const ScalarField& turbulentViscosity() const noexcept
    {
        return turbulenceModel_->turbulentViscosity();
    }

    /// Get wall distance field
    [[nodiscard]] const ScalarField& wallDistance() const noexcept
    {
        return turbulenceModel_->wallDistance();
    }

    /// Get y+ field
    [[nodiscard]] const FaceData<Scalar>& yPlus() const noexcept
    {
        return turbulenceModel_->yPlus();
    }

    /// Get wall shear stress field
    [[nodiscard]] const FaceData<Scalar>& wallShearStress() const noexcept
    {
        return turbulenceModel_->wallShearStress();
    }

    /// Whether the meshWave wall-distance loop converged during initialization
    [[nodiscard]] bool wallDistanceConverged() const noexcept
    {
        return turbulenceModel_ && turbulenceModel_->wallDistanceConverged();
    }

private:

// Private members

// Mesh & Numerical
    const Mesh& mesh_;
    const BoundaryConditions& bcManager_;
    const GradientScheme& gradientScheme_;
    const ConvectionScheme& convectionScheme_;

// Physical properties

    /// Fluid density
    Scalar rho_ = S(1.225);

    /// Dynamic viscosity
    Scalar mu_ = S(1.7894e-5);

    /// Kinematic viscosity
    Scalar nu_ = S(1.7894e-5/1.225);

// Algorithm parameters

    /// Under-relaxation factor for velocity
    Scalar alphaU_ = S(0.7);

    /// Under-relaxation factor for pressure
    Scalar alphaP_= S(0.3);

    /// Under-relaxation factor for k
    Scalar alphaK_ = S(0.5);

    /// Under-relaxation factor for omega
    Scalar alphaOmega_ = S(0.5);

    /// Maximum number of iterations
    int maxIterations_ = 500;

    /// Convergence tolerance
    Scalar tolerance_ = S(1e-3);

/// Turbulence and constraints

    /// Turbulence model
    std::unique_ptr<kOmegaSST> turbulenceModel_ = nullptr;

    /// Field constraint system
    std::unique_ptr<Constraint> constraintSystem_ = nullptr;

// Solution fields

    /// Velocity fields
    ScalarField Ux_;
    ScalarField Uy_;
    ScalarField Uz_;

    /// Pressure field
    ScalarField p_;

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

    /// Pressure-correction gradient
    VectorField gradPCorrPrecomputed_;

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

// Turbulence residual fields

    /// Previous-iteration turbulence fields used for residual computation
    ScalarField kPrev_;
    ScalarField omegaPrev_;

// Matrix and solver fields

    /// Matrix constructor and solver object
    std::unique_ptr<Matrix> matrixConstruct_ = nullptr;

    /// Linear solver for momentum equations
    std::unique_ptr<LinearSolver> momentumSolver_;

    /// Linear solver for pressure correction equation
    std::unique_ptr<LinearSolver> pressureSolver_;

// Residual tracking for convergence

    /// Track pressure correction RMS before reset
    Scalar lastPressureCorrectionRMS_ = S(1e9);

    /// Track turbulence field changes between iterations
    Scalar lastKResidual_ = S(1e9);
    Scalar lastOmegaResidual_ = S(1e9);

    /// First-iteration reference values for scaled residuals
    Scalar massImbalance0_ = S(0.0);
    Scalar velocityResidual0_ = S(0.0);
    Scalar pressureResidual0_ = S(0.0);
    Scalar kResidual0_ = S(0.0);
    Scalar omegaResidual0_ = S(0.0);

// Enable verbose console output

    bool debug_ = false;

// Algorithm steps (called only by solve())

    /// Solve momentum equations for each velocity component
    void solveMomentumEquations();

    /// Update face mass fluxes using Rhie-Chow interpolation
    void updateRhieChowFlowRate();

    /// Assemble and solve the pressure correction Poisson equation
    void solvePressureCorrection();

    /// Apply SIMPLE velocity correction: U = U* - D*gradPCorr
    void correctVelocity();

    /// Update pressure with under-relaxation and reset pCorr
    void correctPressure();

    /// Update face mass fluxes from pressure correction gradient
    void correctFlowRate();

    /// Advance k-omega SST turbulence equations (if enabled)
    void solveTurbulence();

    /// Check convergence against scaled residual tolerance
    [[nodiscard]] bool checkConvergence();

// Private methods

    /// Compute limited velocity gradients into gradUx_/gradUy_/gradUz_
    /// and assemble gradU_
    void updateVelocityGradients();

    /// Compute mass imbalance across domain
    Scalar massImbalance() const noexcept;

    /// Compute velocity residual
    Scalar velocityResidual() const noexcept;

    /// Compute pressure residual
    Scalar pressureResidual() const noexcept;

    /// Add Σf (νEff)f · (∇U)f^T · Sf to momentum source terms
    void addTransposeGradientSource();

    /// Solve a single momentum component equation
    void solveMomentumEquation
    (
        TransportEquation& eq,
        const ScalarField& componentPrev
    );
};
