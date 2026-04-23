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

#include "Scalar.hpp"
#include "Mesh.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"
#include "BoundaryConditions.hpp"
#include "GradientScheme.hpp"
#include "ConvectionSchemes.hpp"
#include "LinearSolvers.hpp"

/// Forward declaration
class Matrix;

class kOmegaSST
{
public:

    /// Constructor
    kOmegaSST
    (
        const Mesh& mesh,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionSchemes& kScheme,
        const ConvectionSchemes& omegaScheme
    );

    /// Copy constructor and assignment
    kOmegaSST(const kOmegaSST&) = delete;
    kOmegaSST& operator=(const kOmegaSST&) = delete;

    /// Move constructor and assignment
    kOmegaSST(kOmegaSST&&) = delete;
    kOmegaSST& operator=(kOmegaSST&&) = delete;

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
     * @param gradU Pre-computed velocity gradient tensor field
     */
    void solve
    (
        const VectorField& U,
        const FaceFluxField& flowRateFace,
        const TensorField& gradU
    );

// Accessor methods

    /**
     * @brief Get turbulent kinetic energy field
     * @return Const reference to k field
     */
    [[nodiscard]] const ScalarField& k() const noexcept { return k_; }

    /**
     * @brief Get specific dissipation rate field
     * @return Const reference to omega field
     */
    [[nodiscard]] const ScalarField& omega() const noexcept { return omega_; }

    /**
     * @brief Get turbulent kinematic viscosity field
     * @return Const reference to nut field
     */
    [[nodiscard]] const ScalarField& turbulentViscosity() const noexcept
    {
        return nut_;
    }

    /**
     * @brief Get wall distance field
     * @return Const reference to wall distance field
     */
    [[nodiscard]] const ScalarField& wallDistance() const noexcept
    {
        return wallDistance_;
    }

    /**
     * @brief Get yPlus field
     * @return Const reference to yPlus field
     */
    [[nodiscard]] const FaceData<Scalar>& yPlus() const noexcept
    {
        return yPlus_;
    }

    /**
     * @brief Get wall shear stress field
     * @return Const reference to wall shear stress field
     */
    [[nodiscard]] const FaceData<Scalar>& wallShearStress() const noexcept
    {
        return wallShearStress_;
    }

    /**
     * @brief Get wall-function nut values on wall faces
     * @return Const reference to nutWallFace field
     */
    [[nodiscard]] const FaceData<Scalar>& nutWall() const noexcept
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
        Scalar a1;           ///< Stress limiter constant
        Scalar c1;           ///< Production limiter coefficient
        Scalar kappa;        ///< Karman constant
        Scalar E;            ///< Wall function constant

    // Numerical constants
        Scalar betaStar;     ///< Destruction coefficient (= Cmu = 0.09)
    };

    static constexpr ModelConstants coeffs_
    {
        .sigmaK1     = S(0.85),
        .sigmaOmega1 = S(0.5),
        .beta1       = S(0.075),
        .gamma1      = S(5.0) / S(9.0),
        .sigmaK2     = S(1.0),
        .sigmaOmega2 = S(0.856),
        .beta2       = S(0.0828),
        .gamma2      = S(0.44),
        .a1          = S(0.31),
        .c1          = S(10.0),
        .kappa       = S(0.41),
        .E           = S(9.8),
        .betaStar    = S(0.09)
    };

private:

// Mesh references

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

    /// Reference to BCs
    const BoundaryConditions& bcManager_;

    /// Reference to gradient scheme
    const GradientScheme& gradientScheme_;

    /// Reference to k convection scheme
    const ConvectionSchemes& kConvectionScheme_;

    /// Reference to omega convection scheme
    const ConvectionSchemes& omegaConvectionScheme_;

// Turbulence fields

    /// Turbulent kinetic energy
    ScalarField k_{S(1e-6)};

    /// Specific dissipation rate
    ScalarField omega_{S(1.0)};

    /// Turbulent kinematic viscosity
    ScalarField nut_;

    /// Distance to nearest wall
    ScalarField wallDistance_{S(1.0)};

    /// Coordinates of the nearest wall point per cell (for mesh-wave propagation)
    VectorField nearestWallPoint_;

    /// Kinematic wall shear stress magnitude (tau/rho) [m^2/s^2]
    FaceData<Scalar> wallShearStress_;

    /// y+
    FaceData<Scalar> yPlus_;


// Auxiliary fields

    /// Wall-function nut values on wall faces (nutkWallFunction)
    FaceData<Scalar> nutWall_;

    /// Dynamic omega wall-function values on faces
    FaceData<Scalar> omegaWallFunctionFaceValues_{
        std::numeric_limits<Scalar>::quiet_NaN()};

    /// Area-based weight per wall face (face area / total wall area of cell)
    FaceData<Scalar> wallFaceWeight_;

    /// Indices into mesh_.faces for faces with wall-function BCs
    std::vector<size_t> wallFunctionFaceIndices_;

    /// Unique cell indices adjacent to wall-function faces
    std::vector<size_t> wallCellIndices_;

    /// Wall-to-total boundary area fraction per wall cell
    std::vector<Scalar> wallCellFraction_;

    /// Area-weighted wall-function omega per wall cell
    std::vector<Scalar> wallCellOmega_;

    /// Precomputed Cmu^0.25 (avoids repeated std::pow calls)
    const Scalar Cmu25_ = std::sqrt(std::sqrt(coeffs_.betaStar));

    /// y+ crossover between viscous sublayer and log region
    Scalar yPlusLam_ = S(11.225);

    /// Minimum wall-normal distance to prevent division by zero
    static constexpr Scalar minWallDist_ = S(1e-20);

// Physical properties

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

// Persistent-state writers

    /// Compute y+ crossover via fixed-point iteration
    void yPlusLam();

    /// Update wall distance field using mesh-wave coordinate propagation
    void updateWallDistance();

    /// Build per-cell corner weights for wall-function faces
    void wallFunctionWeights();

    /// Update dynamic omega wall-function boundary values per face
    void updateOmegaWallFunctionBoundaryValues();

    /// Pre-set wall-cell omega to area-weighted wall-function values
    void applyOmegaWallCellValues();

    /// Update wall-function nut on wall faces (nutkWallFunction)
    void updateNutWallFace();

    /**
     * @brief Update SST turbulent viscosity with limiter
     * @param f23 F2 or F2*F3 blending selector
     * @param strainRateField Strain rate magnitude
     */
    void updateTurbulentViscosity
    (
        const ScalarField& f23,
        const ScalarField& strainRateField
    );

    /**
     * @brief Update wall shear stress using parallel velocity only
     *
     * @details
     * Computes kinematic wall shear stress (tau/rho) [m^2/s^2]:
     * - Viscous sublayer (y+ < yPlusLam): tau/rho = nu * U_tan / y
     * - Log layer (y+ >= yPlusLam): tau/rho = uTau^2 = Cmu^0.5 * k
     *
     * @param U Current velocity field
     */
    void updateWallShearStress(const VectorField& U);

    /// Update y+ field on wall-function faces
    void updateYPlus();

// Transient-quantity producers (return by value)

    /**
     * @brief Compute strain rate magnitude from velocity gradient tensor
     * @param gradU Velocity gradient tensor field
     * @return Strain rate magnitude per cell
     */
    ScalarField strainRate(const TensorField& gradU) const;

    /// Compute limited cell gradient of k
    VectorField gradientK() const;

    /// Compute limited cell gradient of omega
    VectorField gradientOmega() const;

    /**
     * @brief Compute cross-diffusion term from k and omega gradients
     * @param gK Cell gradient of k
     * @param gO Cell gradient of omega
     * @return Cross-diffusion field
     */
    ScalarField crossDiffusion
    (
        const VectorField& gK,
        const VectorField& gO
    ) const;

    /**
     * @brief Compute F1 blending function
     * @param CDkOmega Cross-diffusion field
     * @return F1 blending function
     */
    ScalarField F1(const ScalarField& CDkOmega) const;

    /// Compute F2 blending function
    ScalarField F2() const;

    /// Compute F3 blending switch function (optional)
    ScalarField F3() const;

    /**
     * @brief Compute F23 blending selector
     * @param f2 F2 blending function
     * @param f3 F3 blending switch
     * @return F23 = F2 (or F2*F3 when useF3_ is true)
     */
    ScalarField F23(const ScalarField& f2, const ScalarField& f3) const;

    /**
     * @brief Compute cell velocity divergence from face mass fluxes
     * @param flowRateFace Face volume flow rates
     * @return Cell-centred divergence of U
     */
    ScalarField divU(const FaceFluxField& flowRateFace) const;

// Transport equation solvers

    /// k and omega production term pair
    struct ProductionTerms
    {
        ScalarField k;
        ScalarField omega;
    };

    /**
     * @brief Compute k and omega production terms
     * @param f1 SST F1 blending function
     * @param f23 F2 or F2*F3 blending selector
     * @param strainRateField Strain rate magnitude
     * @return Cell-centred production fields for k and omega
     */
    ProductionTerms productionTerms
    (
        const ScalarField& f1,
        const ScalarField& f23,
        const ScalarField& strainRateField
    ) const;

    /**
     * @brief Override k production at wall-adjacent cells
     * @param U Current velocity field
     * @param productionK k-production field (mutated in place)
     */
    void overrideWallCellProduction
    (
        const VectorField& U,
        ScalarField& productionK
    );

    /**
     * @brief Solve the omega transport equation
     * @param flowRateFace Face volume flow rates
     * @param f1 SST F1 blending function
     * @param productionOmega Omega production field
     * @param CDkOmega Cross-diffusion field
     * @param divUField Cell velocity divergence
     * @param gradOmega Limited cell gradient of omega
     */
    void solveOmegaEquation
    (
        const FaceFluxField& flowRateFace,
        const ScalarField& f1,
        const ScalarField& productionOmega,
        const ScalarField& CDkOmega,
        const ScalarField& divUField,
        const VectorField& gradOmega
    );

    /**
     * @brief Solve the k transport equation
     * @param flowRateFace Face volume flow rates
     * @param f1 SST F1 blending function
     * @param productionK k-production field
     * @param divUField Cell velocity divergence
     * @param gradK Limited cell gradient of k
     */
    void solveKEquation
    (
        const FaceFluxField& flowRateFace,
        const ScalarField& f1,
        const ScalarField& productionK,
        const ScalarField& divUField,
        const VectorField& gradK
    );

// Utilities

    /// Blend two constants using SST blending function
    static Scalar blend(Scalar f, Scalar cInner, Scalar cOuter) noexcept
    {
        return f * (cInner - cOuter) + cOuter;
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
