/******************************************************************************
 * @file KOmegaSST.h
 * @brief k-omega SST (Shear Stress Transport) Turbulence Model Implementation
 * 
 * @details This header file defines the KOmegaSST class, which implements the 
 * k-omega SST turbulence model developed by Menter. The model combines the 
 * robustness of the k-omega model in the near-wall region with the freestream 
 * independence of the k-epsilon model in the far field.
 * 
 * The model solves two transport equations:
 * - Turbulent kinetic energy (k)
 * - Specific dissipation rate (omega)
 * 
 * Key features:
 * - Automatic switching between k-omega and k-epsilon formulations
 * - Enhanced near-wall treatment with proper wall functions
 * - Wall distance calculation using Poisson equation
 * - Proper wall shear stress computation
 * - Production limiting for numerical stability
 * - Support for both steady-state and transient simulations
 * 
 * @see Menter, F.R. (1994). "Two-equation eddy-viscosity turbulence models 
 *       for engineering applications." AIAA Journal, 32(8), 1598-1605.
 */

#ifndef KOMEGASST_H
#define KOMEGASST_H

#include <vector>
#include <memory>

#include "Scalar.h"
#include "Vector.h"
#include "Cell.h"
#include "Face.h"
#include "CellData.h"
#include "FaceData.h"
#include "BoundaryConditions.h"
#include "GradientScheme.h"
#include "LinearSolvers.h"

// Forward declaration
class Matrix;

/**
 * @brief k-omega SST (Shear Stress Transport) Turbulence Model
 * 
 * @details This class implements the k-omega SST turbulence model, which is a 
 * two-equation eddy-viscosity turbulence model that automatically switches 
 * between k-omega and k-epsilon formulations based on local flow conditions.
 * 
 * The SST model is particularly well-suited for:
 * - Adverse pressure gradient flows
 * - Separated flows
 * - Wall-bounded flows
 * - Aerospace applications
 */
class KOmegaSST 
{
public:
    /**
     * @brief Constructor for KOmegaSST turbulence model
     * 
     * @param faces Reference to vector of mesh faces
     * @param cells Reference to vector of mesh cells
     * @param bc Reference to boundary conditions manager
     * @param gradScheme Reference to gradient computation scheme
     * 
     * @details Initializes the turbulence model with mesh data and numerical schemes.
     * Allocates memory for turbulence fields and sets default model constants.
     */
    KOmegaSST
    (
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme
    );
    
    /// Destructor
    ~KOmegaSST();

    /**
     * @brief Initialize turbulence fields with initial conditions
     * 
     * @param U_field Initial velocity field
     * @param nu_lam Laminar kinematic viscosity
     * 
     * @details Sets initial values for k, omega, and other turbulence fields.
     * Uses empirical correlations for initial turbulence levels based on 
     * velocity field and laminar viscosity.
     */
    void initialize(const VectorField& U_field, Scalar nu_lam);

    /**
     * @brief Solve turbulence equations for current iteration
     * 
     * @param U_field Current velocity field
     * @param nu_lam Laminar kinematic viscosity
     * 
     * @details This is the main solve method that:
     * 1. Calculates velocity gradients per cell internally
     * 2. Calculates blending functions F1 and F2
     * 3. Computes production and cross-diffusion terms
     * 4. Solves the k-equation
     * 5. Solves the omega-equation
     * 6. Updates turbulent viscosity
     * 7. Applies wall corrections
     * 
     * @note Velocity gradients are computed internally for efficiency
     */
    void solve
    (
        const VectorField& U_field,
        Scalar nu_lam
    );

    /**
     * @brief Calculate wall distance field using Poisson equation
     * 
     * Solves the Poisson equation ∇²y = 1 with boundary conditions
     * y = 0 at walls to compute the distance to the nearest wall for each cell.
     * This is used for wall function calculations and y+ determination.
     */
    void calculateWallDistance();

    /**
     * @brief Solve the omega transport equation
     * 
     * @param U_field Current velocity field
     * @param gradU Vector of velocity gradient fields
     * @param nu_lam Laminar kinematic viscosity
     * 
     * @details Assembles and solves the linear system for the omega equation
     * using the finite volume method with appropriate boundary conditions.
     */
    void solveOmegaEquation
    (
        const VectorField& U_field,
        const std::vector<VectorField>& gradU,
        Scalar nu_lam
    );

    /**
     * @brief Apply near-wall treatment for omega equation
     * 
     * @param nu_lam Laminar kinematic viscosity
     * 
     * @details Applies wall functions and ensures proper behavior of omega
     * in the viscous sublayer and buffer region.
     */
    void applyNearWallTreatmentOmega(Scalar nu_lam);

    /**
     * @brief Solve the k transport equation
     * 
     * @param U_field Current velocity field
     * @param gradU Vector of velocity gradient fields
     * @param nu_lam Laminar kinematic viscosity
     * 
     * @details Assembles and solves the linear system for the k equation
     * using the finite volume method with appropriate boundary conditions.
     */
    void solveKEquation
    (
        const VectorField& U_field,
        const std::vector<VectorField>& gradU,
        Scalar nu_lam
    );

    /**
     * @brief Calculate turbulent viscosity field
     * 
     * @param U_field Current velocity field
     * @param gradU Vector of velocity gradient fields
     * @param nu_lam Laminar kinematic viscosity
     * 
     * @details Computes μₜ = ρk/ω with appropriate limiting and wall corrections.
     * The turbulent viscosity is limited to prevent unrealistic values and ensure
     * numerical stability.
     */
    void calculateTurbulentViscosity
    (
        const VectorField& U_field,
        const std::vector<VectorField>& gradU,
        Scalar nu_lam
    );

    /**
     * @brief Apply wall corrections for turbulent viscosity
     * 
     * @details Ensures μₜ = 0 at walls and proper behavior in the viscous sublayer.
     * Applies damping functions and ensures smooth transition between wall and 
     * bulk flow regions.
     */
    void applyWallCorrections();

    /**
     * @brief Calculate wall shear stress using parallel velocity only
     * 
     * @param U_field Current velocity field
     * @param nu_lam Laminar kinematic viscosity
     * 
     * @details Computes wall shear stress as:
     * τ_wall = μ_eff * (∂U_parallel/∂n)_wall
     * 
     * This is used for wall function calculations and friction coefficient
     * determination.
     */
    void calculateWallShearStress(const VectorField& U_field, Scalar nu_lam);

    // ==================== GETTER METHODS ====================
    
    /**
     * @brief Get turbulent kinetic energy field
     * @return Const reference to k field
     */
    const ScalarField& getK() const { return k; }
    
    /**
     * @brief Get specific dissipation rate field
     * @return Const reference to omega field
     */
    const ScalarField& getOmega() const { return omega; }
    
    /**
     * @brief Get turbulent viscosity field
     * @return Const reference to μₜ field
     */
    const ScalarField& getTurbulentViscosity() const { return mu_t; }
    
    /**
     * @brief Get wall distance field
     * @return Const reference to wall distance field
     */
    const ScalarField& getWallDistance() const { return wallDistance; }
    
    /**
     * @brief Get wall shear stress field
     * @return Const reference to wall shear stress field
     */
    const ScalarField& getWallShearStress() const { return wallShearStress; }

    /**
     * @brief Get effective viscosity (laminar + turbulent)
     * 
     * @param nu_lam Laminar kinematic viscosity
     * @return Effective viscosity field ν_eff = ν_lam + ν_t
     * 
     * @details Returns the total effective viscosity field that combines
     * laminar and turbulent contributions for use in momentum equations.
     */
    ScalarField getEffectiveViscosity(Scalar nu_lam) const;

    // ==================== SETTER METHODS ====================
    
    /**
     * @brief Set model constants for customization
     * 
     * @param C_mu Turbulent viscosity constant (default: 0.09)
     * @param beta_1 Destruction coefficient for inner region (default: 0.075)
     * @param beta_2 Destruction coefficient for outer region (default: 0.0828)
     * 
     * @details Allows modification of key model constants for specific
     * applications or calibration to experimental data.
     */
    void setModelConstants
    (
        Scalar C_mu = 0.09,
        Scalar beta_1 = 0.075,
        Scalar beta_2 = 0.0828
    );
    
    /**
     * @brief Configure transient simulation parameters
     * 
     * @param enable Whether to enable transient mode
     * @param dt Time step size for transient calculations
     * @param theta Crank-Nicolson parameter (0.5 = implicit, 1.0 = explicit)
     * 
     * @details Sets up the model for transient simulations with appropriate
     * time integration parameters.
     */
    void setTransientMode(bool enable, Scalar dt = 0.001, Scalar theta = 0.5);
    
    /**
     * @brief Update old values for time stepping
     * 
     * @details Called before solve() to store current field values for
     * use in transient calculations. This method should be called at
     * the beginning of each time step.
     */
    void updateOldValues();

private:
    // ==================== TURBULENCE FIELDS ====================
    
    ScalarField k;                  ///< Turbulent kinetic energy [m²/s²]
    ScalarField omega;              ///< Specific dissipation rate [1/s]
    ScalarField mu_t;               ///< Turbulent viscosity [kg/(m·s)]
    ScalarField wallDistance;       ///< Distance to nearest wall [m]
    ScalarField wallShearStress;    ///< Wall shear stress magnitude [Pa]
    
    // ==================== GRADIENT FIELDS ====================
    
    VectorField grad_k;             ///< Gradient of turbulent kinetic energy
    VectorField grad_omega;         ///< Gradient of specific dissipation rate

    // ==================== MESH REFERENCES ====================
    
    const std::vector<Face>& allFaces;      ///< Reference to all mesh faces
    const std::vector<Cell>& allCells;      ///< Reference to all mesh cells
    const BoundaryConditions& bcManager;    ///< Reference to boundary conditions
    const GradientScheme& gradientScheme;   ///< Reference to gradient scheme

    // ==================== AUXILIARY FIELDS ====================
    
    ScalarField F1;                 ///< Blending function 1 (k-omega to k-epsilon)
    ScalarField F2;                 ///< Blending function 2 (viscous sublayer)
    ScalarField production_k;       ///< Production of turbulent kinetic energy
    ScalarField production_omega;   ///< Production of specific dissipation rate
    ScalarField cross_diffusion;    ///< Cross-diffusion term in omega equation

    // ==================== MODEL CONSTANTS ====================
    
    /**
     * @brief Structure containing all model constants
     * 
     * @details The SST model uses different constants for the inner (k-omega)
     * and outer (k-epsilon) regions, with blending functions to smoothly
     * transition between them.
     */
    struct ModelConstants {
        // k-omega model constants (inner region)
        Scalar sigma_k1 = 0.85;     ///< k diffusion coefficient (inner)
        Scalar sigma_omega1 = 0.5;  ///< omega diffusion coefficient (inner)
        Scalar beta_1 = 0.075;      ///< Destruction coefficient (inner)
        Scalar gamma_1 = 5.0/9.0;   ///< Production coefficient (inner)
        
        // k-epsilon model constants (outer region)
        Scalar sigma_k2 = 1.0;      ///< k diffusion coefficient (outer)
        Scalar sigma_omega2 = 0.856; ///< omega diffusion coefficient (outer)
        Scalar beta_2 = 0.0828;     ///< Destruction coefficient (outer)
        Scalar gamma_2 = 0.44;      ///< Production coefficient (outer)
        
        // Common constants
        Scalar C_mu = 0.09;         ///< Turbulent viscosity constant
        Scalar a1 = 0.31;           ///< Stress limiter constant
        Scalar kappa = 0.41;        ///< Karman constant
        Scalar E = 9.8;             ///< Wall function constant
        
        // Numerical constants
        Scalar beta_star = 0.09;    ///< Modified destruction coefficient
        Scalar sigma_d = 2.0 * sigma_omega2; ///< Cross-diffusion coefficient
    } constants;

    // ==================== TRANSIENT PARAMETERS ====================
    
    bool enableTransient;        ///< Flag to enable transient simulation
    Scalar dt;                  ///< Time step size [s]
    Scalar theta;               ///< Crank-Nicolson parameter (0.5 = implicit)
    
    // ==================== PREVIOUS TIME STEP FIELDS ====================
    
    ScalarField k_old;          ///< k from previous time step
    ScalarField omega_old;      ///< omega from previous time step
    
    // ==================== PHYSICAL PROPERTIES ====================
    
    Scalar rho;                 ///< Fluid density [kg/m³]

    // ==================== NUMERICAL TOOLS ====================
    
    std::unique_ptr<Matrix> matrixConstruct;  ///< Matrix constructor for equation solving

    // ==================== PRIVATE HELPER METHODS ====================
    
    /**
     * @brief Calculate blending functions F1 and F2
     * 
     * @param gradU Vector of velocity gradient fields
     * @param nu_lam Laminar kinematic viscosity
     * 
     * @details Computes the blending functions that determine the transition
     * between k-omega and k-epsilon formulations based on local flow conditions.
     */
    void calculateBlendingFunctions
    (
        const std::vector<VectorField>& gradU,
        Scalar nu_lam
    );

    /**
     * @brief Calculate production terms for k and omega equations
     * 
     * @param gradU Vector of velocity gradient fields
     * @param rho Fluid density
     * 
     * @details Computes P_k = μₜ(∇U:∇U) and P_ω = γ/νₜ * P_k with
     * appropriate limiting to prevent unrealistic values.
     */
    void calculateProductionTerms(const std::vector<VectorField>& gradU);

    /**
     * @brief Calculate cross-diffusion term for omega equation
     * 
     * @param rho Fluid density
     * 
     * @details Computes the cross-diffusion term that appears only in the
     * k-epsilon formulation and is blended out in the k-omega region.
     */
    void calculateCrossDiffusion();

    /**
     * @brief Apply turbulence boundary conditions
     * 
     * @param fieldName Name of the turbulence field ("k" or "omega")
     * @param field Reference to the field to apply BCs to
     * @param U_field Current velocity field
     * @param nu_lam Laminar kinematic viscosity
     * 
     * @details Applies appropriate boundary conditions for turbulence fields:
     * - Wall: k = 0, omega from wall functions
     * - Inlet: k and omega from turbulence intensity
     * - Outlet: Zero gradient
     */
    void applyTurbulenceBoundaryConditions
    (
        const std::string& fieldName,
        ScalarField& field,
        const VectorField& U_field,
        Scalar nu_lam
    );

    /**
     * @brief Calculate y+ value for wall treatment
     * 
     * @param cellIdx Cell index
     * @param U_field Velocity field
     * @param nu_lam Laminar kinematic viscosity
     * @param rho Density
     * @return y+ value (dimensionless wall distance)
     * 
     * @details Computes y+ = y*u_τ/ν where y is wall distance, u_τ is
     * friction velocity, and ν is kinematic viscosity. Used to determine
     * appropriate wall function treatment.
     */
    Scalar calculateYPlus
    (
        size_t cellIdx,
        const VectorField& U_field,
        Scalar nu_lam
    ) const;

    /**
     * @brief Limit production to prevent unrealistic values
     * 
     * @param P_k Production of k
     * @param rho Density
     * @param cellIdx Cell index
     * @return Limited production value
     * 
     * @details Applies production limiting to prevent k from growing
     * unrealistically large, which can cause numerical instability.
     * The limiting is based on the local turbulence timescale.
     */
    Scalar limitProduction(Scalar P_k, size_t cellIdx) const;
};

#endif // KOMEGASST_H