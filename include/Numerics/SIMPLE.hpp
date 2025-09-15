/******************************************************************************
 * @file SIMPLE.h
 * @brief SIMPLE algorithm for incompressible Navier-Stokes equations
 * 
 * This class implements the Semi-Implicit Method for Pressure-Linked 
 * Equations (SIMPLE) for solving incompressible flow on unstructured finite
 * volume meshes. The algorithm handles velocity-pressure coupling through 
 * pressure correction and includes support for k-omega SST turbulence 
 * modeling with wall functions.
 * 
 * @class SIMPLE
 * 
 * The SIMPLE class provides a complete CFD solver featuring:
 * - Pressure-velocity coupling via SIMPLE algorithm
 * - Rhie-Chow interpolation for collocated grid arrangement  
 * - Momentum equations with implicit under-relaxation
 * - Pressure correction with mass conservation enforcement
 * - k-omega SST turbulence model with wall distance calculation
 * - Field constraints for numerical stability
 * 
 * Key algorithmic features:
 * - Segregated solution of momentum and pressure equations
 * - Face velocity reconstruction preventing checkerboard oscillations
 * - Deferred correction for higher-order convection schemes
 * - Adaptive convergence monitoring and residual tracking
 *****************************************************************************/

#ifndef SIMPLE_H
#define SIMPLE_H

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
#include "KOmegaSST.hpp"
#include "Constraint.hpp"


/**
 * @brief SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) solver
 * 
 * This class implements the SIMPLE algorithm for solving incompressible
 * Navier-Stokes equations using finite volume discretization. It handles
 * velocity-pressure coupling through pressure correction and includes
 * support for turbulence modeling.
 */
class SIMPLE 
{
public:
    /**
     * @brief Constructor for SIMPLE solver
     * @param faces Reference to face data
     * @param cells Reference to cell data
     * @param bc Reference to boundary conditions
     * @param gradScheme Reference to gradient scheme
     * @param convScheme Reference to convection scheme
     */
    SIMPLE
    (
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells, 
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionScheme& convScheme
    );

    /**
     * @brief Main SIMPLE algorithm execution
     */
    void solve();
    
    /**
     * @brief Initialize solution fields
     * @param initialVelocity Initial velocity field
     * @param initialPressure Initial pressure field
     */
    void initialize(const Vector& initialVelocity, Scalar initialPressure);
    
    /**
     * @brief Solve momentum equations for velocity components
     */
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
    
    /**
     * @brief Correct pressure field using pressure correction
     */
    void correctPressure();
    
    /**
     * @brief Correct mass fluxes using updated velocity field
     */
    void correctFlowRate();
    
    /**
     * @brief Solve turbulence equations if turbulence modeling is enabled
     * 
     * Solves the turbulent transport equations (k-omega SST) to update
     * turbulent viscosity field. Only executes if turbulence is enabled.
     */
    void solveTurbulence();
    
    /**
     * @brief Check if solution has converged
     * @return true if converged, false otherwise
     */
    bool checkConvergence();

    /**
     * @brief Get velocity field
     * @return Reference to velocity field
     */
    const VectorField& getVelocity() const { return U; }
    
    /**
     * @brief Get pressure field
     * @return Reference to pressure field
     */
    const ScalarField& getPressure() const { return p; }
    
    /**
     * @brief Get Rhie-Chow face flow rate field
     * @return Reference to Rhie-Chow flow rate field
     */
    const FaceFluxField& getRhieChowFlowRate() const { return RhieChowFlowRate; }
    
    /**
     * @brief Get x-component of velocity field
     * @return Scalar field of x-velocity components
     */
    ScalarField getVelocityX() const;
    
    /**
     * @brief Get y-component of velocity field
     * @return Scalar field of y-velocity components
     */
    ScalarField getVelocityY() const;
    
    /**
     * @brief Get z-component of velocity field
     * @return Scalar field of z-velocity components
     */
    ScalarField getVelocityZ() const;

    /**
     * @brief Set under-relaxation factors
     * @param alpha_U Velocity under-relaxation factor
     * @param alpha_p Pressure under-relaxation factor
     */
    void setRelaxationFactors(Scalar alpha_U, Scalar alpha_p);
    
    /**
     * @brief Set convergence tolerance
     * @param tol Convergence tolerance value
     */
    void setConvergenceTolerance(Scalar tol);
    
    /**
     * @brief Set maximum number of iterations
     * @param maxIter Maximum iteration count
     */
    void setMaxIterations(int maxIter);
    
    /**
     * @brief Enable or disable turbulence modeling
     * @param enable Whether to enable turbulence modeling
     */
    void enableTurbulenceModeling(bool enable = true);
    
    /**
     * @brief Set physical properties
     * @param rho_new Fluid density
     * @param mu_new Dynamic viscosity
     */
    void setPhysicalProperties(Scalar rho_new, Scalar mu_new) 
    {
        this->rho = rho_new;
        this->mu = mu_new;
        this->nu = mu_new / rho_new;
    }
    
    /**
     * @brief Get constraint system pointer
     * @return Pointer to constraint system
     */
    Constraint* getConstraintSystem() { return constraintSystem.get(); }
    
    /**
     * @brief Get turbulent kinetic energy field
     * @return Pointer to k field (null if turbulence disabled)
     */
    const ScalarField* getTurbulentKineticEnergy() const;
    
    /**
     * @brief Get specific dissipation rate field
     * @return Pointer to omega field (null if turbulence disabled)
     */
    const ScalarField* getSpecificDissipationRate() const;
    
    /**
     * @brief Get turbulent viscosity field
     * @return Pointer to nu_t field (null if turbulence disabled)
     */
    const ScalarField* getTurbulentViscosity() const;
    
    /**
     * @brief Get wall distance field
     * @return Pointer to wall distance field (null if turbulence disabled)
     */
    const ScalarField* getWallDistance() const;

private:

    /// Mesh references
    const std::vector<Face>& allFaces;
    const std::vector<Cell>& allCells;
    const BoundaryConditions& bcManager;
    const GradientScheme& gradientScheme;
    const ConvectionScheme& convectionScheme;

    /// Physical properties
    Scalar rho;             ///< Fluid density
    Scalar mu;              ///< Dynamic viscosity
    Scalar nu;              ///< Kinematic viscosity

    /// Algorithm parameters
    Scalar alpha_U;         ///< Under-relaxation factor for velocity
    Scalar alpha_p;         ///< Under-relaxation factor for pressure
    int maxIterations;      ///< Maximum number of iterations
    Scalar tolerance;       ///< Convergence tolerance
    bool enableTurbulence;  ///< Enable turbulence modeling
    
    
    /// Turbulence model
    std::unique_ptr<KOmegaSST> turbulenceModel;
    
    /// Field constraint system
    std::unique_ptr<Constraint> constraintSystem;

    /// Velocity field
    VectorField U;
    
    /// Pressure field
    ScalarField p;
    
    /// Pressure correction field
    ScalarField pCorr;
    
    /// Track pressure correction RMS before reset
    Scalar lastPressureCorrectionRMS = S(1e9);

    /// Velocity from previous iteration
    VectorField U_prev;
    
    /// X-component from previous iteration
    ScalarField U_x_prev;
    
    /// Y-component from previous iteration
    ScalarField U_y_prev;
    
    /// Z-component from previous iteration
    ScalarField U_z_prev;

    /// Face velocity field (current iteration)
    FaceVectorField U_f_avg;
    
    /// Face velocity from previous iteration
    FaceVectorField U_f_avg_prev;
    
    /// Mass flux through faces (Rhie-Chow)
    FaceFluxField RhieChowFlowRate;
    
    /// Mass flux from previous iteration
    FaceFluxField RhieChowFlowRate_prev;

    /// Cell diffusion coefficients for momentum
    ScalarField DU;
    
    /// Face diffusion coefficients for momentum
    FaceFluxField DUf;

    // Gradient fields
    /// Pressure gradient field
    VectorField gradP;
    
    /// Pressure correction gradient field
    VectorField gradPCorr;

    /// Matrix constructor and solver object
    std::unique_ptr<Matrix> matrixConstruct;
    
    /**
     * @brief Calculate mass imbalance across domain
     * @return Mass imbalance value
     */
    Scalar calculateMassImbalance() const;
    
    /**
     * @brief Calculate velocity residual
     * @return Velocity residual value
     */
    Scalar calculateVelocityResidual() const;
    
    /**
     * @brief Calculate pressure residual
     * @return Pressure residual value
     */
    Scalar calculatePressureResidual() const;    
};

#endif