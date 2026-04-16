/******************************************************************************
 * @file SIMPLE.hpp
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

#include <array>
#include <vector>
#include <memory>

#include "Mesh.hpp"
#include "BoundaryConditions.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"
#include "GradientScheme.hpp"
#include "ConvectionScheme.hpp"
#include "Matrix.hpp"
#include "LinearSolvers.hpp"
#include "ErrorHandler.hpp"
#include "kOmegaSST.hpp"
#include "Constraint.hpp"


class SIMPLE
{
public:

    /**
     * @brief Constructor for SIMPLE solver
     * @param mesh Mesh view (nodes, faces, cells)
     * @param bc Reference to boundary conditions
     * @param gradScheme Reference to gradient scheme
     * @param convSchemes Reference to per-equation convection schemes
     */
    SIMPLE
    (
        const Mesh& mesh,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionSchemes& convSchemes
    );

    /// Destructor (defined for forward-declared types)
    ~SIMPLE();

    /// Main SIMPLE algorithm execution
    void solve();

    /**
     * @brief Initialize solution fields
     * @param initialVelocity Initial velocity field
     * @param initialPressure Initial pressure field
     * @param initialK Initial turbulent kinetic energy
     * @param initialOmega Initial specific dissipation rate
     * @param enableTurbulence Enable k-omega SST turbulence modeling
     */
    void initialize
    (
        const Vector& initialVelocity,
        Scalar initialPressure,
        Scalar initialK = S(1e-6),
        Scalar initialOmega = S(1.0),
        bool enableTurbulence = false
    );

// Accessor methods

    /**
     * @brief Get velocity field
     * @return Const reference to velocity vector field
     */
    const VectorField& velocity() const noexcept { return U_; }

    /**
     * @brief Get pressure field
     * @return Const reference to pressure scalar field
     */
    const ScalarField& pressure() const noexcept { return p_; }

// Setter methods

    /**
     * @brief Set under-relaxation factors
     * @param alphaU Velocity under-relaxation factor
     * @param alphaP Pressure under-relaxation factor
     * @param alphaK k under-relaxation factor (default: 0.5)
     * @param alphaOmega omega under-relaxation factor (default: 0.5)
     */
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

    /**
     * @brief Set convergence tolerance
     * @param tol Convergence tolerance for residuals
     */
    void setConvergenceTolerance(Scalar tol) noexcept { tolerance_ = tol; }

    /**
     * @brief Set maximum number of iterations
     * @param maxIter Maximum number of SIMPLE iterations
     */
    void setMaxIterations(int maxIter) noexcept { maxIterations_ = maxIter; }

    /**
     * @brief Enable or disable verbose console output
     * @param d True to enable debug output
     */
    void setDebug(bool d) noexcept { debug_ = d; }

    /**
     * @brief Set physical properties
     * @param rho Fluid density
     * @param mu Dynamic viscosity
     */
    void setPhysicalProperties(Scalar rho, Scalar mu)
    {
        if (rho <= S(0.0))
        {
            FatalError("setPhysicalProperties: density must be positive");
        }
        rho_ = rho;
        mu_  = mu;
        nu_  = mu / rho;
    }

    /**
     * @brief Set linear solver for momentum equations
     * @param solver Configured LinearSolver for momentum
     */
    void setMomentumSolver(const LinearSolver& solver)
    {
        momentumSolver_ = solver;
    }

    /**
     * @brief Set linear solver for pressure correction equation
     * @param solver Configured LinearSolver for pressure
     */
    void setPressureSolver(const LinearSolver& solver)
    {
        pressureSolver_ = solver;
    }

    /**
     * @brief Set linear solvers for turbulence equations
     * @param kSolver Configured LinearSolver for k equation
     * @param omegaSolver Configured LinearSolver for omega equation
     */
    void setTurbulenceSolvers
    (
        const LinearSolver& kSolver,
        const LinearSolver& omegaSolver
    );

// Getter methods

    /**
     * @brief Get constraint system
     * @return Reference to constraint system
     */
    [[nodiscard("Constraint system needed for applying fixed-value cells")]]
    Constraint& constraintSystem() noexcept { return *constraintSystem_; }

    /**
     * @brief Get turbulent kinetic energy field
     * @return Reference to k field
     */
    const ScalarField& turbulentKineticEnergy() const noexcept
    {
        return turbulenceModel_->k();
    }

    /**
     * @brief Get specific dissipation rate field
     * @return Reference to omega field
     */
    const ScalarField& specificDissipationRate() const noexcept
    {
        return turbulenceModel_->omega();
    }

    /**
     * @brief Get turbulent viscosity field
     * @return Reference to nut field
     */
    const ScalarField& turbulentViscosity() const noexcept
    {
        return turbulenceModel_->turbulentViscosity();
    }

    /**
     * @brief Get wall distance field
     * @return Reference to wall distance field
     */
    const ScalarField& wallDistance() const noexcept
    {
        return turbulenceModel_->wallDistance();
    }

    /**
     * @brief Get y+ field
     * @return Reference to y+ field
     */
    const FaceData<Scalar>& yPlus() const noexcept
    {
        return turbulenceModel_->yPlus();
    }

    /**
     * @brief Get wall shear stress field
     * @return Reference to wall shear stress field
     */
    const FaceData<Scalar>& wallShearStress() const noexcept
    {
        return turbulenceModel_->wallShearStress();
    }

private:

// Private members

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;
    const BoundaryConditions& bcManager_;
    const GradientScheme& gradientScheme_;
    const ConvectionSchemes& convectionSchemes_;

// Physical properties

    /// Fluid density
    Scalar rho_ = S(1.225);
    /// Dynamic viscosity
    Scalar mu_ = S(1.7894e-5);
    /// Kinematic viscosity
    Scalar nu_ = S(1.7894e-5/1.225);

/// Algorithm parameters

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

    /// Enable verbose console output
    bool debug_ = false;

    /// Turbulence model
    std::unique_ptr<kOmegaSST> turbulenceModel_ = nullptr;

    /// Field constraint system
    std::unique_ptr<Constraint> constraintSystem_ = nullptr;

    /// Velocity field
    VectorField U_;

    /// Pressure field
    ScalarField p_;

    /// Pressure correction field
    ScalarField pCorr_;

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

    /// Velocity from previous iteration
    VectorField UPrev_;

    /// Face velocity field (current iteration)
    FaceVectorField UAvgf_;

    /// Face velocity from previous iteration
    FaceVectorField UAvgPrevf_;

    /// Mass flux through faces (Rhie-Chow)
    FaceFluxField RhieChowFlowRate_;

    /// Mass flux from previous iteration
    FaceFluxField RhieChowFlowRatePrev_;

    /// Cell diffusion coefficients for momentum
    ScalarField DU_;

    /// Face diffusion coefficients for momentum
    FaceFluxField DUf_;

// Gradient fields

    /// Pressure gradient field
    VectorField gradP_;

    /// Pressure correction gradient field
    VectorField gradPCorr_;

    /// Velocity gradients (shared between momentum,
    /// transpose source, turbulence)
    std::array<VectorField, 3> gradU_;

    /// Matrix constructor and solver object
    std::unique_ptr<Matrix> matrixConstruct_ = nullptr;

    /// Linear solver for momentum equations
    LinearSolver momentumSolver_ = LinearSolver("momentum", S(1e-6), 1000);

    /// Linear solver for pressure correction equation
    LinearSolver pressureSolver_ = LinearSolver("pCorr", S(1e-6), 1000);

// Private types

    /// Velocity component selector used in place of raw char/int indices
    enum class Component : int { X = 0, Y = 1, Z = 2 };

// Algorithm steps (called only by solve())

    /// Solve momentum equations for each velocity component
    void solveMomentumEquations();

    /// Calculate face mass fluxes using Rhie-Chow interpolation
    void calculateRhieChowFlowRate();

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

    /**
     * @brief Extract a single component from a VectorField into a ScalarField
     * @param V Source vector field
     * @param component Component to extract (X, Y, or Z)
     * @return ScalarField containing the extracted component
     */
    static ScalarField extractComponent
    (
        const VectorField& V,
        Component component
    );

    /// Calculate mass imbalance across domain
    Scalar calculateMassImbalance() const;

    /// Calculate velocity residual
    Scalar calculateVelocityResidual() const;

    /// Calculate pressure residual
    Scalar calculatePressureResidual() const;

    /**
     * @brief Calculate transpose gradient source term for momentum equations
     *
     * @details Computes the explicit source term: ∇·(ν_eff · (∇U)^T)
     * This term arises from the full viscous stress tensor τ = μ(∇U + (∇U)^T)
     * and is non-zero when viscosity varies spatially (turbulent flows).
     *
     * Implemented as face-based divergence:
     * Σ_f (ν_eff)_f · (∇U)_f^T · S_f
     *
     * @param nuEffFace Face-based effective viscosity (ν + ν_t)
     * @param transposeSourceX Output: x-momentum source term
     * @param transposeSourceY Output: y-momentum source term
     * @param transposeSourceZ Output: z-momentum source term
     */
    void calculateTransposeGradientSource
    (
        const FaceData<Scalar>& nuEffFace,
        ScalarField& transposeSourceX,
        ScalarField& transposeSourceY,
        ScalarField& transposeSourceZ
    );

    /**
     * @brief Solve a single momentum component equation
     *
     * @details Builds and solves the discretized momentum equation for one
     * velocity component. Handles matrix assembly, under-relaxation, and
     * linear solver iteration.
     *
     * @param component Component to solve (X, Y, or Z)
     * @param eq Transport equation data for this component
     * @param componentPrev Previous iteration velocity component
     */
    void solveMomentumComponent
    (
        Component component,
        TransportEquation& eq,
        const ScalarField& componentPrev
    );

    /**
     * @brief Build gradient transpose vector from velocity gradients
     *
     * @details Constructs the transpose gradient columns:
     * - X: [∂u/∂x, ∂v/∂x, ∂w/∂x]
     * - Y: [∂u/∂y, ∂v/∂y, ∂w/∂y]
     * - Z: [∂u/∂z, ∂v/∂z, ∂w/∂z]
     *
     * @param gradUx Gradient of x-velocity
     * @param gradUy Gradient of y-velocity
     * @param gradUz Gradient of z-velocity
     * @param component Which column of the transpose to build
     * @return Vector representing the transpose gradient column
     */
    Vector buildGradientTransposeColumn
    (
        const Vector& gradUx,
        const Vector& gradUy,
        const Vector& gradUz,
        Component component
    ) const noexcept;
};
