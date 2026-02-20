/******************************************************************************
 * @file SIMPLE.hpp
 * @brief SIMPLE algorithm for incompressible Navier-Stokes equations
 * 
 * This header implements the Semi-Implicit Method for Pressure-Linked 
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

#ifndef SIMPLE_HPP
#define SIMPLE_HPP

#include <iostream>
#include <vector>
#include <memory>

#include "Face.hpp"
#include "Cell.hpp"
#include "BoundaryConditions.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"
#include "GradientScheme.hpp"
#include "ConvectionScheme.hpp"
#include "Matrix.hpp"
#include "LinearSolvers.hpp"
#include "kOmegaSST.hpp"
#include "Constraint.hpp"


class SIMPLE 
{
public:

    /**
     * @brief Constructor for SIMPLE solver
     * @param faces Reference to face data
     * @param cells Reference to cell data
     * @param bc Reference to boundary conditions
     * @param gradScheme Reference to gradient scheme
     * @param convSchemes Reference to per-equation convection schemes
     */
    SIMPLE
    (
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionSchemes& convSchemes
    );

    /// Main SIMPLE algorithm execution
    void solve();

    /**
     * @brief Initialize solution fields
     * @param initialVelocity Initial velocity field
     * @param initialPressure Initial pressure field
     * @param initialK Initial turbulent kinetic energy
     * @param initialOmega Initial specific dissipation rate
     */
    void initialize
    (
        const Vector& initialVelocity,
        Scalar initialPressure,
        Scalar initialK = S(1e-6),
        Scalar initialOmega = S(1.0)
    );

    /// Solve momentum equations for velocity components
    void solveMomentumEquations();

    /**
     * @brief Calculate mass fluxes using Rhie-Chow interpolation
     *
     * Implements the Rhie-Chow interpolation method to prevent checkerboard
     * pressure oscillations in collocated grids.
     *
     * The face velocity is computed as:
     *
     * U_f = U_f_average + DUf * (∇p_face_average - ∇p_face_accurate)
     *       + (1 - alpha_U) * (U_f_prev - U_f_prev_linear)
     *
     * where D_f is the face diffusion coefficient calculated using the
     * momentum equation diagonals.
     */
    void calculateRhieChowFlowRate();

    /**
     * @brief Solve pressure correction equation
     *
     * Assembles and solves the pressure correction equation:
     * ∇·(DUf ∇p') = -∇·(ρU*)
     */
    void solvePressureCorrection();

    /**
     * @brief Correct velocity field using pressure correction
     *
     * Applies velocity correction as:
     * U_corrected = U* - DU * ∇p'
     *
     * where DU = diag(DUx, DUy, DUz) are the cell-centered diffusion
     * coefficients computed from momentum matrix diagonals: D = V/a_momentum
     */
    void correctVelocity();

    /// Correct pressure field using pressure correction
    void correctPressure();

    /// Correct mass fluxes using updated velocity field
    void correctFlowRate();

    /**
     * @brief Solve turbulence equations if turbulence modeling is enabled
     */
    void solveTurbulence();

    /// Check if solution has converged
    bool checkConvergence();

// Accessor methods

    /// Get velocity field
    const VectorField& getVelocity() const { return U_; }

    /// Get pressure field
    const ScalarField& getPressure() const { return p_; }

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
        Scalar alphaK = 0.5,
        Scalar alphaOmega = 0.5
    )
    {
        alphaU_ = alphaU;
        alphaP_ = alphaP;
        alphaK_ = alphaK;
        alphaOmega_ = alphaOmega;
    }

    /// Set convergence tolerance
    void setConvergenceTolerance(Scalar tol)
    {
        tolerance_ = tol;
    }

    /// Set maximum number of iterations
    void setMaxIterations(int maxIter)
    {
        maxIterations_ = maxIter;
    }

    /// Enable or disable turbulence modeling
    void enableTurbulenceModeling(bool enable = true)
    {
        enableTurbulence_ = enable;

        if (enable)
        {
            std::cout
                << "k-omega SST turbulence modeling enabled."
                << std::endl;
        }
        else
        {
            std::cout
                << "Laminar flow (turbulence modeling"
                << " disabled)."
                << std::endl;
        }
    }

    /**
     * @brief Set physical properties
     * @param rhoNew Fluid density
     * @param muNew Dynamic viscosity
     */
    void setPhysicalProperties(Scalar rhoNew, Scalar muNew)
    {
        rho_ = rhoNew;
        mu_ = muNew;
        nu_ = muNew / rhoNew;
    }

    /// Set linear solver for momentum equations
    void setMomentumSolver(const LinearSolver& solver)
    {
        momentumSolver_ = solver;
    }

    /// Set linear solver for pressure correction equation
    void setPressureSolver(const LinearSolver& solver)
    {
        pressureSolver_ = solver;
    }

    /**
     * @brief Set linear solvers for turbulence equations
     * @param kSolver Configured LinearSolver for k equation
     * @param omegaSolver Configured LinearSolver for omega equation
     *
     * Must be called after initialize() when turbulence is enabled.
     */
    void setTurbulenceSolvers
    (
        const LinearSolver& kSolver,
        const LinearSolver& omegaSolver
    )
    {
        if (turbulenceModel_)
        {
            turbulenceModel_->kSolverSettings() = kSolver;
            turbulenceModel_->omegaSolverSettings()
                = omegaSolver;
        }
    }

    /// Get constraint system pointer
    Constraint* getConstraintSystem() { return constraintSystem_.get(); }

    /// Get turbulent kinetic energy field (null if turbulence disabled)
    const ScalarField* getTurbulentKineticEnergy() const
    {
        if (enableTurbulence_ && turbulenceModel_)
        {
            return &(turbulenceModel_->getk());
        }
        return nullptr;
    }

    /// Get specific dissipation rate field (null if turbulence disabled)
    const ScalarField* getSpecificDissipationRate() const
    {
        if (enableTurbulence_ && turbulenceModel_)
        {
            return &(turbulenceModel_->getOmega());
        }
        return nullptr;
    }

    /// Get turbulent viscosity field (null if turbulence disabled)
    const ScalarField* getTurbulentViscosity() const
    {
        if (enableTurbulence_ && turbulenceModel_)
        {
            return &(turbulenceModel_
                ->getTurbulentViscosity());
        }
        return nullptr;
    }

    /// Get wall distance field (null if turbulence disabled)
    const ScalarField* getWallDistance() const
    {
        if (enableTurbulence_ && turbulenceModel_)
        {
            return &(turbulenceModel_
                ->getWallDistance());
        }
        return nullptr;
    }

private:

// Private members

    /// Mesh references
    const std::vector<Face>& allFaces_;
    const std::vector<Cell>& allCells_;
    const BoundaryConditions& bcManager_;
    const GradientScheme& gradientScheme_;
    const ConvectionSchemes& convectionScheme_;

    /// Physical properties

    /// Fluid density
    Scalar rho_;
    /// Dynamic viscosity
    Scalar mu_;
    /// Kinematic viscosity
    Scalar nu_;

    /// Algorithm parameters

    /// Under-relaxation factor for velocity
    Scalar alphaU_;
    /// Under-relaxation factor for pressure
    Scalar alphaP_;
    /// Under-relaxation factor for k
    Scalar alphaK_;
    /// Under-relaxation factor for omega
    Scalar alphaOmega_;
    /// Maximum number of iterations
    int maxIterations_;
    /// Convergence tolerance
    Scalar tolerance_;
    /// Enable turbulence modeling
    bool enableTurbulence_;

    /// Turbulence model
    std::unique_ptr<kOmegaSST> turbulenceModel_;

    /// Field constraint system
    std::unique_ptr<Constraint> constraintSystem_;

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
    std::vector<VectorField> gradU_;

    /// Matrix constructor and solver object
    std::unique_ptr<Matrix> matrixConstruct_;

    /// Linear solver for momentum equations
    LinearSolver momentumSolver_;

    /// Linear solver for pressure correction equation
    LinearSolver pressureSolver_;

// Private methods

    /**
     * @brief Extract a component (x=0, y=1, z=2) from a VectorField
     * @param V Source vector field
     * @param component Component index (0=x, 1=y, 2=z)
     * @param name Name for the resulting scalar field
     * @return ScalarField containing the extracted component
     */
    static ScalarField extractComponent
    (
        const VectorField& V,
        int component,
        const std::string& name
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
     * Computes the explicit source term: ∇·(ν_eff · (∇U)^T)
     * This term arises from the full viscous stress tensor τ = μ(∇U + (∇U)^T)
     * and is non-zero when viscosity varies spatially (turbulent flows).
     *
     * Implemented as face-based divergence following OpenFOAM's approach:
     * Σ_f (ν_eff)_f · (∇U)_f^T · S_f
     *
     * @param nu_eff Effective viscosity field (ν + ν_t)
     * @param gradUx Pre-computed gradient of x-velocity
     * @param gradUy Pre-computed gradient of y-velocity
     * @param gradUz Pre-computed gradient of z-velocity
     * @param transposeSourceX Output: x-momentum source term
     * @param transposeSourceY Output: y-momentum source term
     * @param transposeSourceZ Output: z-momentum source term
     */
    void calculateTransposeGradientSource
    (
        const ScalarField& nu_eff,
        const VectorField& gradUx,
        const VectorField& gradUy,
        const VectorField& gradUz,
        ScalarField& transposeSourceX,
        ScalarField& transposeSourceY,
        ScalarField& transposeSourceZ
    );

    /**
     * @brief Solve a single momentum component equation
     *
     * Builds and solves the discretized momentum equation for one velocity
     * component (x, y, or z). Handles matrix assembly, under-relaxation,
     * and linear solver iteration.
     *
     * @param component Component identifier ('x', 'y', or 'z')
     * @param U_component Current velocity component field
     * @param source Source term for this component (includes pressure gradient)
     * @param U_component_prev Previous iteration velocity component
     * @param nu_eff Effective viscosity field
     * @param grad_component Pre-computed cell gradients
     *        of this velocity component
     */
    void solveMomentumComponent
    (
        char component,
        ScalarField& U_component,
        const ScalarField& source,
        const ScalarField& U_component_prev,
        const ScalarField& nu_eff,
        const VectorField& grad_component
    );

    /**
     * @brief Build gradient transpose vector from velocity gradients
     *
     * Constructs the transpose gradient columns for a given component:
     * - component 0 (x): [∂u/∂x, ∂v/∂x, ∂w/∂x]
     * - component 1 (y): [∂u/∂y, ∂v/∂y, ∂w/∂y]
     * - component 2 (z): [∂u/∂z, ∂v/∂z, ∂w/∂z]
     *
     * @param gradUx Gradient of x-velocity
     * @param gradUy Gradient of y-velocity
     * @param gradUz Gradient of z-velocity
     * @param component Which column of transpose (0=x, 1=y, 2=z)
     * @return Vector representing the transpose gradient column
     */
    Vector buildGradientTransposeColumn
    (
        const Vector& gradUx,
        const Vector& gradUy,
        const Vector& gradUz,
        int component
    ) const;
};

#endif // SIMPLE_HPP