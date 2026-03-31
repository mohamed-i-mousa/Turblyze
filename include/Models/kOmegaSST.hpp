/******************************************************************************
 * @file kOmegaSST.hpp
 * @brief k-omega SST Shear Stress Transport Turbulence Model Implementation
 *
 * @details This header file defines the KOmegaSST class, which implements the
 * k-omega SST turbulence model developed by Menter. The model combines the
 * robustness of the k-omega model in the near-wall region with the freestream
 * independence of the k-epsilon model in the far field.
 *
 * @class kOmegaSST
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

#pragma once

#include <vector>
#include <memory>
#include <span>

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

/// Forward declaration
class Matrix;

class kOmegaSST
{
public:

    /// Constructor
    kOmegaSST
    (
        std::span<const Face> faces,
        std::span<const Cell> cells,
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
     *
     * @details
     * This is the main solve method that:
     * 1. Computes k and omega gradients
     * 2. Calculates blending functions F1 and F2
     * 3. Computes strain rate, production, and cross-diffusion terms
     * 4. Solves the omega-equation with under-relaxation
     * 5. Solves the k-equation with under-relaxation
     * 6. Updates turbulent viscosity and wall shear stress
     *
     * @param U Current velocity field
     * @param flowRateFace Face volume flow rates
     * @param gradU Pre-computed velocity gradients {gradUx, gradUy, gradUz}
     */
    void solve
    (
        const VectorField& U,
        const FaceFluxField& flowRateFace,
        std::span<const VectorField> gradU
    );

// Accessor methods

    /**
     * @brief Get turbulent kinetic energy field
     * @return Const reference to k field
     */
    [[nodiscard("Turbulence field needed for output or boundary treatment")]]
    const ScalarField& k() const noexcept { return k_; }

    /**
     * @brief Get specific dissipation rate field
     * @return Const reference to omega field
     */
    [[nodiscard("Turbulence field needed for output or boundary treatment")]]
    const ScalarField& omega() const noexcept
    {
        return omega_;
    }

    /**
     * @brief Get turbulent kinematic viscosity field
     * @return Const reference to nut field
     */
    [[nodiscard("Turbulence field needed for output or boundary treatment")]]
    const ScalarField& turbulentViscosity() const noexcept
    {
        return nut_;
    }

    /**
     * @brief Get wall distance field
     * @return Const reference to wall distance field
     */
    [[nodiscard("Turbulence field needed for output or boundary treatment")]]
    const ScalarField& wallDistance() const noexcept
    {
        return wallDistance_;
    }

    /**
     * @brief Get yPlus field
     * @return Const reference to yPlus field
     */
    [[nodiscard("Turbulence field needed for output or boundary treatment")]]
    const FaceData<Scalar>& yPlus() const noexcept
    {
        return yPlus_;
    }

    /**
     * @brief Get wall shear stress field
     * @return Const reference to wall shear stress field
     */
    [[nodiscard("Turbulence field needed for output or boundary treatment")]]
    const FaceData<Scalar>& wallShearStress() const noexcept
    {
        return wallShearStress_;
    }

    /**
     * @brief Get effective viscosity (laminar + turbulent)
     * @return Effective viscosity field nu_eff = nu_lam + nu_t
     */
    [[nodiscard("Computed field required for momentum equation diffusion")]]
    ScalarField effectiveViscosity() const;

    /**
     * @brief Get wall-function nut values on wall faces
     * @return Const reference to nutWallFace field
     */
    [[nodiscard("Turbulence field needed for output or boundary treatment")]]
    const FaceData<Scalar>& nutWall() const noexcept
    {
        return nutWall_;
    }

    /**
     * @brief Access k-equation linear solver for external configuration
     * @return Reference to k linear solver
     */
    LinearSolver& kSolverSettings() noexcept { return kSolver_; }

    /**
     * @brief Access omega-equation linear solver for external configuration
     * @return Reference to omega linear solver
     */
    LinearSolver& omegaSolverSettings() noexcept { return omegaSolver_; }

    /**
     * @brief Enable or disable verbose console output
     * @param d True to enable debug output
     */
    void setDebug(bool d) noexcept { debug_ = d; }

// Model constants

    /// Structure containing all model constants
    struct ModelConstants
    {
    // k-omega model constants (near-wall region)
        Scalar sigmaK1;      ///< k diffusion coefficient
        Scalar sigmaOmega1;  ///< omega diffusion coefficient
        Scalar beta1;        ///< Destruction coefficient
        Scalar gamma1;       ///< Production coefficient

    // k-epsilon model constants (outer region)
        Scalar sigmaK2;      ///< k diffusion coefficient
        Scalar sigmaOmega2;  ///< omega diffusion coefficient
        Scalar beta2;        ///< Destruction coefficient
        Scalar gamma2;       ///< Production coefficient

    // Common constants
        Scalar Cmu;          ///< Turbulent viscosity constant
        Scalar a1;           ///< Stress limiter constant
        Scalar b1;           ///< F2 multiplier in nut denominator
        Scalar c1;           ///< Production limiter coefficient
        Scalar kappa;        ///< Karman constant
        Scalar E;            ///< Wall function constant

    // Numerical constants
        Scalar betaStar;     ///< Modified destruction coefficient
    };

    static constexpr ModelConstants const_
    {
        .sigmaK1     = S(0.85),
        .sigmaOmega1 = S(0.5),
        .beta1       = S(0.075),
        .gamma1      = S(5.0) / S(9.0),
        .sigmaK2     = S(1.0),
        .sigmaOmega2 = S(0.856),
        .beta2       = S(0.0828),
        .gamma2      = S(0.44),
        .Cmu         = S(0.09),
        .a1          = S(0.31),
        .b1          = S(1.0),
        .c1          = S(10.0),
        .kappa       = S(0.41),
        .E           = S(9.8),
        .betaStar    = S(0.09)
    };

private:

// Mesh references

    /// Reference to all mesh faces
    std::span<const Face> allFaces_;

    /// Reference to all mesh cells
    std::span<const Cell> allCells_;

    /// Reference to BCs
    const BoundaryConditions& bcManager_;

    /// Reference to gradient scheme
    const GradientScheme& gradientScheme_;

    /// Reference to k convection scheme
    const ConvectionScheme& kConvectionScheme_;

    /// Reference to omega convection scheme
    const ConvectionScheme& omegaConvectionScheme_;

// Turbulence fields

    /// Turbulent kinetic energy
    ScalarField k_ = ScalarField("k", allCells_.size(), S(1e-6));

    /// Specific dissipation rate
    ScalarField omega_ = ScalarField("omega", allCells_.size(), S(1.0));

    /// Turbulent kinematic viscosity
    ScalarField nut_ = ScalarField("nut", allCells_.size(), S(0.0));

    /// Distance to nearest wall
    ScalarField wallDistance_ =
        ScalarField("wallDistance", allCells_.size(), S(1.0));

    /// Wall shear stress magnitude
    FaceData<Scalar> wallShearStress_ =
        FaceData<Scalar>("wallShearStress", allFaces_.size(), S(0.0));

    /// y+
    FaceData<Scalar> yPlus_ =
        FaceData<Scalar>("yPlus", allFaces_.size(), S(0.0));

// Gradient fields

    /// Gradient of k
    VectorField gradK_ = VectorField("gradK", allCells_.size(), Vector{});

    /// Gradient of Omega
    VectorField gradOmega_ =
        VectorField("gradOmega", allCells_.size(), Vector{});

// Auxiliary fields

    /// Cached velocity divergence (computed once per solve() call)
    ScalarField divU_ = ScalarField("divU", allCells_.size(), S(0.0));

    /// Blending function 1 (k-ω to k-ε)
    ScalarField F1_ = ScalarField("F1", allCells_.size(), S(0.0));

    /// Blending function 2 (viscous sublayer)
    ScalarField F2_ = ScalarField("F2", allCells_.size(), S(0.0));

    /// Blending function 23
    ScalarField F23_ = ScalarField("F23", allCells_.size(), S(0.0));

    /// Blending function 3
    ScalarField F3_ = ScalarField("F3", allCells_.size(), S(0.0));

    /// Production of k
    ScalarField productionK_ =
        ScalarField("productionK", allCells_.size(), S(0.0));

    /// Production of omega
    ScalarField productionOmega_ =
        ScalarField("productionOmega", allCells_.size(), S(0.0));

    /// Strain rate magnitude (shared between production and viscosity)
    ScalarField strainRate_ =
        ScalarField("strainRate", allCells_.size(), S(0.0));

    /// Cross-diffusion term: 2·σω2·(∇k·∇ω)/ω
    ScalarField CDkOmega_ = ScalarField("CDkOmega", allCells_.size(), S(0.0));

    /// Wall-function nut values on wall faces (nutkWallFunction)
    FaceData<Scalar> nutWall_ =
        FaceData<Scalar>("nutWall", allFaces_.size(), S(0.0));

    /// Dynamic omega wall-function values on faces
    FaceData<Scalar> omegaWallFunctionFaceValues_ =
        FaceData<Scalar>
        (
            "omegaWallFunctionFaceValues",
            allFaces_.size(),
            std::numeric_limits<Scalar>::quiet_NaN()
        );

    /// Area-based weight per wall face (face area / total wall area of cell)
    FaceData<Scalar> wallFaceWeight_ =
        FaceData<Scalar>("wallFaceWeight", allFaces_.size(), S(0.0));

    /// Indices into allFaces_ for faces with wall-function BCs
    std::vector<size_t> wallFunctionFaceIndices_;

    /// Unique cell indices adjacent to wall-function faces
    std::vector<size_t> wallCellIndices_;

    /// Wall-to-total boundary area fraction per wall cell
    std::vector<Scalar> wallCellFraction_;

    /// Area-weighted wall-function omega per wall cell
    std::vector<Scalar> wallCellOmega_;

    /// Precomputed Cmu^0.25 (avoids repeated std::pow calls)
    const Scalar Cmu25_ = std::sqrt(std::sqrt(const_.Cmu));

// Physical properties

    /// Fluid density (used for wall shear stress)
    Scalar rho_ = S(1.225);

    /// Laminar kinematic viscosity
    Scalar nu_ = S(1.7894e-5 / 1.225);

    /// Optional SST F3 switch
    bool useF3_ = false;

    /// Positive floor for k
    Scalar kMin_ = S(smallValue);

    /// Positive floor for omega
    Scalar omegaMin_ = S(smallValue);

    /// Maximum turbulent-to-laminar viscosity ratio for omega bound
    Scalar nutMaxCoeff_ = S(1e5);

    /// y+ crossover between viscous and log layers
    Scalar yPlusLam_ = S(11.225);

    /// Under-relaxation factor for k equation
    Scalar alphaK_ = S(0.5);

    /// Under-relaxation factor for omega equation
    Scalar alphaOmega_ = S(0.5);

    /// Enable verbose console output
    bool debug_ = false;

// Numerical tools

    /// Matrix constructor
    std::unique_ptr<Matrix> matrixConstruct_;

    /// Linear solver for k equation
    LinearSolver kSolver_ = LinearSolver("k", S(1e-8), 1000);

    /// Linear solver for omega equation
    LinearSolver omegaSolver_ = LinearSolver("omega", S(1e-8), 1000);

// Private helper methods

    /// Compute y+ crossover via fixed-point iteration
    void calcYPlusLam();

    /// Calculate wall distance field using (Dijkstra's shortest path)
    void calculateWallDistance();

    /// Build per-cell corner weights for wall-function faces
    void buildWallFunctionWeights();

    /**
     * @brief Solve the omega transport equation
     * @param flowRateFace Face volume flow rates
     */
    void solveOmegaEquation(const FaceFluxField& flowRateFace);

    /// Update dynamic omega wall-function boundary values per face
    void updateOmegaWallFunctionBoundaryValues();

    /// Pre-set wall-cell omega to area-weighted wall-function values
    void applyOmegaWallCellValues();

    /**
     * @brief Solve the k transport equation
     * @param flowRateFace Face volume flow rates
     */
    void solveKEquation(const FaceFluxField& flowRateFace);

    /// Update wall-function nut on wall faces (nutkWallFunction)
    void updateNutWallFace();

    /// Calculate SST turbulent viscosity with limiter
    void calculateTurbulentViscosity();

    /**
     * @brief Calculate wall shear stress using parallel velocity only
     *
     * @details
     * Computes wall shear stress:
     * - Viscous sublayer (y+ < yPlusLam): tau = rho * nu * U_tan / y
     * - Log layer (y+ >= yPlusLam): tau = rho * uTau^2 = rho * Cmu^0.5 * k
     *
     * @param U Current velocity field
     */
    void calculateWallShearStress(const VectorField& U);

    /// Calculate blending functions F1 and F2
    void calculateBlendingFunctions();

    /// Calculate production terms for k and omega equations
    void calculateProductionTerms();

    /**
     * @brief Override k production at wall-adjacent cells
     * @param U Current velocity field
     */
    void overrideWallCellProduction(const VectorField& U);

    /**
     * @brief Calculate strain rate magnitude field from velocity gradients
     * @param gradU Vector of velocity gradient fields
     */
    void calculateStrainRate(std::span<const VectorField> gradU);

    /**
     * @brief Calculate gradients of k and/or omega fields
     * @param computeGradK Whether to compute k gradients
     * @param computeGradOmega Whether to compute omega gradients
     */
    void calculateGradients(bool computeGradK, bool computeGradOmega);

    /// Compute cross-diffusion term CDkOmega_ from k and omega gradients
    void calculateCrossDiffusion();

    /// Compute cell velocity divergence from face mass fluxes into divU_
    void computeDivU(const FaceFluxField& flowRateFace);

    /// F3 blending switch function (optional)
    void F3();

    /// F23 blending selector
    void F23();

    /// Calculate y+ value for wall treatment
    void calculateYPlus();

    /// Blend two constants using SST blending function
    Scalar blend(Scalar F1, Scalar c1, Scalar c2) const noexcept
    {
        return F1 * (c1 - c2) + c2;
    }

    /// Bound k and omega to physically meaningful positive values
    void boundTurbulenceFields();

    /// Bound omega with viscosity-ratio limit
    void boundOmega();

    /// Bound k to minimum value
    void boundK();

    /// Log min/max/mean for k, omega, nut fields
    void logFieldDiagnostics() const;
};
