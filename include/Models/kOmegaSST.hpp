/******************************************************************************
 * @file kOmegaSST.hpp
 * @brief k-omega SST Shear Stress Transport Turbulence Model Implementation
 * 
 * This header file defines the KOmegaSST class, which implements the 
 * k-omega SST turbulence model developed by Menter. The model combines the 
 * robustness of the k-omega model in the near-wall region with the freestream
 * independence of the k-epsilon model in the far field.
 * 
 * The model solves two transport equations:
 * - Turbulent kinetic energy
 * - Specific dissipation rate
 * 
 * Key features:
 * - Automatic switching between k-omega and k-epsilon formulations
 * - Enhanced near-wall treatment with proper wall functions
 * - Wall distance calculation using iterative propagation approach
 * - Proper wall shear stress computation
 * - Production limiting for numerical stability
 *****************************************************************************/


#ifndef K_OMEGA_SST_HPP
#define K_OMEGA_SST_HPP

#include <vector>
#include <memory>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Cell.hpp"
#include "Face.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"
#include "BoundaryConditions.hpp"
#include "GradientScheme.hpp"
#include "ConvectionScheme.hpp"
#include "LinearSolvers.hpp"

// Forward declaration
class Matrix;

class kOmegaSST 
{
public:

    /// Constructor
    kOmegaSST
    (
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionScheme& kScheme,
        const ConvectionScheme& omegaScheme
    );
    
    /// Destructor
    ~kOmegaSST();

    /**
     * @brief Initialize turbulence fields with initial conditions
     * @param nu Laminar kinematic viscosity
     * @param initialK Initial turbulent kinetic energy
     * @param initialOmega Initial specific dissipation rate
     * @param alphaK Under-relaxation factor for k equation
     * @param alphaOmega Under-relaxation factor for omega equation
     */
    void initialize
    (
        Scalar nu,
        Scalar initialK,
        Scalar initialOmega,
        Scalar alphaK,
        Scalar alphaOmega
    );

    /**
     * @brief Solve turbulence equations for current iteration
     * @param U Current velocity field
     * @param flowRateFace Face volume flow rates
     * @param gradU Pre-computed velocity gradients {gradUx, gradUy, gradUz}
     *
     * @details This is the main solve method that:
     * 1. Computes k and omega gradients
     * 2. Calculates blending functions F1 and F2
     * 3. Computes strain rate, production, and cross-diffusion terms
     * 4. Solves the omega-equation with under-relaxation
     * 5. Solves the k-equation with under-relaxation
     * 6. Updates turbulent viscosity and wall shear stress
     */
    void solve
    (
        const VectorField& U,
        const FaceFluxField& flowRateFace,
        const std::vector<VectorField>& gradU
    );

// Accessor methods

    /**
     * @brief Get turbulent kinetic energy field
     * @return Const reference to k field
     */
    const ScalarField& getk() const { return k_; }

    /**
     * @brief Get specific dissipation rate field
     * @return Const reference to omega field
     */
    const ScalarField& getOmega() const { return omega_; }

    /**
     * @brief Get turbulent kinematic viscosity field
     * @return Const reference to νₜ field
     */
    const ScalarField& getTurbulentViscosity() const { return nut_; }

    /**
     * @brief Get wall distance field
     * @return Const reference to wall distance field
     */
    const ScalarField& getWallDistance() const { return wallDistance_; }

    /**
     * @brief Get wall shear stress field
     * @return Const reference to wall shear stress field
     */
    const ScalarField& getWallShearStress() const { return wallShearStress_; }

    /**
     * @brief Get effective viscosity (laminar + turbulent)
     * @return Effective viscosity field ν_eff = ν_lam + ν_t
     */
    ScalarField getEffectiveViscosity() const;

    /// Access k-equation linear solver for external configuration
    LinearSolver& kSolverSettings() { return kSolver_; }

    /// Access omega-equation linear solver for external configuration
    LinearSolver& omegaSolverSettings() { return omegaSolver_; }

    /**
     * @brief Enable or disable verbose console output
     * @param d True to enable debug output
     */
    void setDebug(bool d) { debug_ = d; }

private:

// Turbulence fields

    /// Turbulent kinetic energy [m²/s²]
    ScalarField k_;

    /// Specific dissipation rate [1/s]
    ScalarField omega_;

    /// Turbulent kinematic viscosity [m²/s]
    ScalarField nut_;

    /// Distance to nearest wall [m]
    ScalarField wallDistance_;

    /// Wall shear stress magnitude [Pa]
    ScalarField wallShearStress_;

// Gradient fields

    /// Gradient of k
    VectorField gradK_;

    /// Gradient of Omega
    VectorField gradOmega_;

// Mesh references

    /// Reference to all mesh faces
    const std::vector<Face>& allFaces_;

    /// Reference to all mesh cells
    const std::vector<Cell>& allCells_;

    /// Reference to BCs
    const BoundaryConditions& bcManager_;

    /// Reference to gradient scheme
    const GradientScheme& gradientScheme_;

    /// Reference to k convection scheme
    const ConvectionScheme& kConvectionScheme_;

    /// Reference to omega convection scheme
    const ConvectionScheme& omegaConvectionScheme_;

// Auxiliary fields

    /// Blending function 1 (k-ω to k-ε)
    ScalarField F1_;

    /// Blending function 2 (viscous sublayer)
    ScalarField F2_;

    /// Production of k
    ScalarField productionK_;

    /// Production of ω
    ScalarField productionOmega_;

    /// Cross-diffusion term in ω equation
    ScalarField crossDiffusion_;

    /// Strain rate magnitude (shared between production and viscosity)
    ScalarField strainRate_;

    /// Dynamic omega wall-function values on faces
    FaceData<Scalar> omegaWallFunctionFaceValues_;

    /// Number of omega-wall-function faces touching each owner cell
    std::vector<size_t> omegaWallFaceCountPerCell_;

    /// Corner weight per cell (1 / omegaWallFaceCountPerCell)
    std::vector<Scalar> omegaWallCellWeight_;

    /// Number of cells that touch more than one omega-wall face
    size_t omegaMultiWallCellCount_ = 0;

// Model constants
    
    /// Structure containing all model constants
    struct ModelConstants
    {
    // k-omega model constants (inner region)

        /// k diffusion coefficient (inner)
        Scalar sigmaK1 = 0.85;

        /// omega diffusion coefficient (inner)
        Scalar sigmaOmega1 = 0.5;

        /// Destruction coefficient (inner)
        Scalar beta1 = 0.075;

        /// Production coefficient (inner)
        Scalar gamma1 = 5.0/9.0;

    // k-epsilon model constants (outer region)

        /// k diffusion coefficient (outer)
        Scalar sigmaK2 = 1.0;

        /// omega diffusion coefficient (outer)
        Scalar sigmaOmega2 = 0.856;

        /// Destruction coefficient (outer)
        Scalar beta2 = 0.0828;

        /// Production coefficient (outer)
        Scalar gamma2 = 0.44;

    // Common constants

        /// Turbulent viscosity constant
        Scalar Cmu = 0.09;

        /// Stress limiter constant
        Scalar a1 = 0.31;

        /// F2 multiplier in nut denominator
        Scalar b1 = 1.0;

        /// Production limiter coefficient
        Scalar c1 = 10.0;
        
        /// Karman constant
        Scalar kappa = 0.41;
        
        /// Wall function constant
        Scalar E = 9.8;

    // Numerical constants
        
        /// Modified destruction coefficient
        Scalar betaStar = 0.09;
        
        /// Cross-diffusion coefficient
        Scalar sigmaD = 2.0 * sigmaOmega2;
    } constants_;

// Physical properties

    /// Fluid density [kg/m³]
    Scalar rho_;

    /// Laminar kinematic viscosity [m²/s]
    Scalar nu_;

    /// Optional SST F3 switch
    bool useF3_;

    /// Positive floor for k
    Scalar kMin_;

    /// Positive floor for omega
    Scalar omegaMin_;

    /// Under-relaxation factor for k equation
    Scalar alphaK_;

    /// Under-relaxation factor for omega equation
    Scalar alphaOmega_;

    /// Enable verbose console output
    bool debug_ = false;

// Numerical tools

    /// Matrix constructor
    std::unique_ptr<Matrix> matrixConstruct_;

    /// Linear solver for k equation
    LinearSolver kSolver_;

    /// Linear solver for omega equation
    LinearSolver omegaSolver_;

// Private helper methods

    /**
     * @brief Calculate wall distance field using an iterative propagation
     * approach (Dijkstra's shortest path)
     */
    void calculateWallDistance();

    /**
     * @brief Build per-cell corner weights for omega wall-function faces
     */
    void buildOmegaWallCellWeights();

    /**
     * @brief Solve the omega transport equation
     * @param flowRateFace Face volume flow rates
     */
    void solveOmegaEquation(const FaceFluxField& flowRateFace);

    /**
     * @brief Update dynamic omega wall-function boundary values per face
     */
    void updateOmegaWallFunctionBoundaryValues();

    /**
     * @brief Solve the k transport equation
     * @param flowRateFace Face volume flow rates
     */
    void solveKEquation(const FaceFluxField& flowRateFace);

    /**
     * @brief Calculate SST turbulent viscosity with limiter
     */
    void calculateTurbulentViscosity();

    /**
     * @brief Calculate wall shear stress using parallel velocity only
     * @param U Current velocity field
     * @details Computes wall shear stress as:
     * τ_wall = μ_eff * (∂U_parallel/∂n)_wall
     */
    void calculateWallShearStress(const VectorField& U);

    /**
     * @brief Calculate blending functions F1 and F2
     */
    void calculateBlendingFunctions();

    /**
     * @brief Calculate production terms for k and omega equations
     *
     * Computes Pk = μₜ(∇U:∇U) and P_ω = γ/νₜ * P_k
     */
    void calculateProductionTerms();

    /**
     * @brief Override k production at wall-adjacent cells
     *   G = (ν + νₜ) Cμ^0.25 √k |U∥| / (κ y²)
     * @param U Current velocity field
     */
    void overrideWallCellProduction(const VectorField& U);

    /**
     * @brief Calculate strain rate magnitude field from velocity gradients
     * @param gradU Vector of velocity gradient fields
     */
    void calculateStrainRate(const std::vector<VectorField>& gradU);

    /**
     * @brief Calculate gradients of k and/or omega fields
     * @param computeGradK Whether to compute k gradients
     * @param computeGradOmega Whether to compute omega gradients
     */
    void calculateGradients(bool computeGradK, bool computeGradOmega);

    /**
     * @brief Calculate cross-diffusion term for omega equation
     */
    void calculateCrossDiffusion();

    /**
     * @brief F3 blending switch function (optional)
     * @param cellIdx Cell index
     * @return F3 value in [0, 1]
     */
    Scalar F3(size_t cellIdx) const;

    /**
     * @brief F23 blending selector
     * @param cellIdx Cell index
     * @return F23 = F2 (default) or F2*F3 if F3 switch is enabled
     */
    Scalar F23(size_t cellIdx) const;

    /**
     * @brief Calculate y+ value for wall treatment
     * @param cellIdx Cell index
     * @return y+ value (dimensionless wall distance)
     */
    Scalar calculateYPlus(size_t cellIdx) const;

    /**
     * @brief Limit production to prevent unrealistic values
     * @param Pk Production of k
     * @param cellIdx Cell index
     * @return Limited production value
     */
    Scalar limitProduction(Scalar Pk, size_t cellIdx) const;

    /// Blend two constants using SST blending function
    inline Scalar blend
    (
        Scalar F1,
        Scalar c1,
        Scalar c2
    ) const
    {
        return F1 * (c1 - c2) + c2;
    }

    /// Bound k and omega to physically meaningful positive values
    void boundTurbulenceFields();

    /// Log min/max/mean for k, omega, nut fields
    void logFieldDiagnostics() const;
};

#endif // K_OMEGA_SST_HPP
