/******************************************************************************
 * @file kOmegaSST.h
 * @brief k-omega SST Shear Stress Transport Turbulence Model Implementation
 *
 * @details This header file defines the KOmegaSST class, which implements the
 * k-omega SST turbulence model developed by Menter. The model combines the
 * robustness of the k-omega model in the near-wall region with the freestream
 * independence of the k-epsilon model in the far field.
 *
 * @class kOmegaSST
 * The model solves two transport equations:
 * - Turbulent kinetic energy
 * - Specific dissipation rate
 *
 * Features:
 * - Automatic switching between k-omega and k-epsilon formulations
 * - Enhanced near-wall treatment with proper wall functions
 * - Wall distance calculation using iterative propagation approach
 * - Proper wall shear stress computation
 * - Production limiting for numerical stability
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

/// Standard library headers
#include <vector>
#include <memory>
#include <limits>

/// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "Mesh.h"
#include "CellData.h"
#include "FaceData.h"
#include "BoundaryConditions.h"
#include "GradientScheme.h"
#include "ConvectionSchemes.h"
#include "LinearSolvers.h"

// *************************** Forward Declarations ***************************

class Matrix;

// ****************************** class kOmegaSST *****************************

class kOmegaSST
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    kOmegaSST
    (
        const Mesh& mesh,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionSchemes& kScheme,
        LinearSolver& kSolver,
        const ConvectionSchemes& omegaScheme,
        LinearSolver& omegaSolver,
        const Scalar nu,
        const Scalar initialK,
        const Scalar initialOmega,
        const Scalar alphaK,
        const Scalar alphaOmega,
        const bool debug
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    kOmegaSST(const kOmegaSST&) = delete;
    kOmegaSST& operator=(const kOmegaSST&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    kOmegaSST(kOmegaSST&&) = delete;
    kOmegaSST& operator=(kOmegaSST&&) = delete;

    /// Destructor
    ~kOmegaSST() noexcept;

// ***************************** Solve kOmegaSST *****************************

    /// Solve turbulence equations for current iteration
    void solve
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz,
        const FaceFluxField& flowRateFace,
        const TensorField& gradU
    );

// ***************************** Accessor Methods *****************************

    /// Get turbulent kinetic energy field
    [[nodiscard]] const ScalarField& k() const noexcept
    {
        return k_;
    }

    /// Get specific dissipation rate field
    [[nodiscard]] const ScalarField& omega() const noexcept
    {
        return omega_;
    }

    /// Get normalised k change from the most recent solve
    [[nodiscard]] Scalar lastKResidual() const noexcept
    {
        return lastKResidual_;
    }

    /// Get normalised omega change from the most recent solve
    [[nodiscard]] Scalar lastOmegaResidual() const noexcept
    {
        return lastOmegaResidual_;
    }

    /// Get turbulent kinematic viscosity field
    [[nodiscard]] const ScalarField& turbulentViscosity() const noexcept
    {
        return nut_;
    }

    /// Get wall distance field
    [[nodiscard]] const ScalarField& wallDistance() const noexcept
    {
        return wallDistance_;
    }

    /// Whether the meshWave wall-distance loop converged
    [[nodiscard]] bool wallDistanceConverged() const noexcept
    {
        return wallDistanceConverged_;
    }

    /// Get yPlus field
    [[nodiscard]] const FaceData<Scalar>& yPlus() const noexcept
    {
        return yPlus_;
    }

    /// Get wall shear stress field
    [[nodiscard]] const FaceData<Scalar>& wallShearStress() const noexcept
    {
        return wallShearStress_;
    }

    /// Get wall-function nut values on wall faces
    [[nodiscard]] const FaceData<Scalar>& nutWall() const noexcept
    {
        return nutWall_;
    }

// ***************************** Model Constants *****************************

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

// *********************** Inlet Condition Calculators ***********************

    /// Calculate inlet/initial turbulent kinetic energy
    [[nodiscard]] static Scalar inletK
    (
        const Vector& velocity,
        Scalar turbulenceIntensity
    ) noexcept;

    /// Calculate inlet/initial specific dissipation rate
    [[nodiscard]] static Scalar inletOmega
    (
        Scalar k,
        Scalar hydraulicDiameter
    ) noexcept;

// ****************************** Private Members *****************************

private:

// Dependencies and services

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

    /// Reference to BCs
    const BoundaryConditions& bcManager_;

    /// Reference to gradient scheme
    const GradientScheme& gradientScheme_;

    /// Reference to k convection scheme
    const ConvectionSchemes& kConvectionScheme_;

    /// Linear solver for k equation
    LinearSolver& kSolver_;

    /// Reference to omega convection scheme
    const ConvectionSchemes& omegaConvectionScheme_;

    /// Linear solver for omega equation
    LinearSolver& omegaSolver_;

    /// Matrix constructor
    std::unique_ptr<Matrix> matrixConstruct_;

// Physical and algorithm parameters

    /// Laminar kinematic viscosity
    Scalar nu_;

    /// Under-relaxation factor for k equation
    Scalar alphaK_;

    /// Under-relaxation factor for omega equation
    Scalar alphaOmega_;

    /// Enable verbose console output
    bool debug_;

    /// Optional SST F3 switch
    bool useF3_ = false;

    /// Maximum turbulent-to-laminar viscosity ratio for omega bound
    Scalar maxViscosityRatio_ = S(1e5);

// Turbulence fields

    /// Turbulent kinetic energy
    ScalarField k_{S(1e-6)};

    /// Specific dissipation rate
    ScalarField omega_{S(1.0)};

    /// Turbulent kinematic viscosity
    ScalarField nut_{S(0.0)};

// Residual tracking

    /// Previous-iteration k/omega snapshots for residual computation
    ScalarField kPrev_;
    ScalarField omegaPrev_;

    /// Normalised k/omega change from the most recent solve
    Scalar lastKResidual_ = S(1e9);
    Scalar lastOmegaResidual_ = S(1e9);

// Production and model-term fields

    /// k production term
    ScalarField Pk_;

    /// Omega production term
    ScalarField POmega_;

    /// Strain-rate magnitude per cell: ||S|| = sqrt(2 S_ij S_ij)
    ScalarField strainRateMag_;

    /// Velocity divergence per cell
    ScalarField divUField_;

    /// Cross-diffusion term per cell
    ScalarField CDkOmega_;

    /// SST blending functions
    ScalarField f1_;
    ScalarField f2_;
    ScalarField f3_;
    ScalarField f23_;

    /// Effective diffusivities for the k and omega transport equations
    ScalarField GammaK_;
    ScalarField GammaOmega_;

    /// Source-term scratch (kept zero; consumed by Matrix::buildMatrix)
    ScalarField kSource_;
    ScalarField omegaSource_;

    /// Cell gradients of k and omega
    VectorField gradK_;
    VectorField gradOmega_;

// Wall distance and wall-function state

    /// Distance to nearest wall
    ScalarField wallDistance_{S(1.0)};

    /// Coordinates of the nearest wall point per cell (for mesh-wave)
    VectorField nearestWallPoint_;

    /// meshWave wall-distance loop convergence flag
    bool wallDistanceConverged_ = false;

    /// Kinematic wall shear stress magnitude (tau/rho) [m^2/s^2]
    FaceData<Scalar> wallShearStress_;

    /// Owner-cell to wall-face perpendicular distance
    FaceData<Scalar> y_;

    /// y+
    FaceData<Scalar> yPlus_;

    /// Wall-function nut values on wall faces (nutkWallFunction)
    FaceData<Scalar> nutWall_;

    /// Dynamic omega wall-function values on faces
    FaceData<Scalar> omegaWall_
    {
        std::numeric_limits<Scalar>::quiet_NaN()
    };

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

    /// Per-cell accumulator of wall-function production contributions
    ScalarField wallProductionAccum_;

    /// Per-cell flag (0/1) marking cells touched by a wall override
    std::vector<char> hasWallOverride_;

    /// Per-cell accumulator of area-weighted wall-function omega values
    std::vector<Scalar> omegaAccum_;

    /// Precomputed Cmu^0.25 (avoids repeated std::pow calls)
    const Scalar Cmu25_ = std::sqrt(std::sqrt(coeffs_.betaStar));

    /// y+ crossover between viscous sublayer and log region
    Scalar yPlusLam_ = S(11.225);

// ***************************** Private Methods *****************************

// Utility helpers

    /// Blend two constants using SST blending function
    static Scalar blend(Scalar f, Scalar cInner, Scalar cOuter) noexcept
    {
        return f * (cInner - cOuter) + cOuter;
    }

// Initialization helpers

    /// Compute y+ crossover via fixed-point iteration
    void yPlusLam();

    /// Update wall distance field using mesh-wave coordinate propagation
    void updateWallDistance();

    /// Build per-cell corner weights for wall-function faces
    void wallFunctionWeights();

    /// Initialize turbulent viscosity with k/omega estimate
    void initializeTurbulentViscosity();

// Wall-function helpers

    /// Update y+ field on wall-function faces
    void updateYPlus();

    /// Update wall-function nut on wall faces (nutkWallFunction)
    void updateNutWall();

    /// Update dynamic omega wall-function boundary values per face
    void updateOmegaWallValues();

    /// Pre-set wall-cell omega to area-weighted wall-function values
    void applyOmegaWallCellValues();

    /// Override k production at wall-adjacent cells (mutates Pk_)
    void overrideWallCellProduction
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz
    );

// Per-iteration model terms

    /// Compute strain rate magnitude into strainRateMag_
    void strainRate(const TensorField& gradU);

    /// Compute cell velocity divergence into divUField_
    void divU(const FaceFluxField& flowRateFace);

    /// Compute k production term into Pk_ (uses strainRateMag_)
    void kProduction();

    /// Compute cross-diffusion term into CDkOmega_ (uses gradK_/gradOmega_)
    void crossDiffusion();

    /// Compute F1 blending function into f1_ (uses CDkOmega_)
    void F1();

    /// Compute F2 blending function into f2_
    void F2();

    /// Compute F3 blending switch function (optional) into f3_
    void F3();

    /// Compute F23 blending selector into f23_ (uses f2_/f3_)
    void F23();

    /// Compute omega production into POmega_ (uses f1_/strainRateMag_)
    void omegaProduction();

    /// Limit Pk_/POmega_ in place using f1_/f23_/strainRateMag_
    void limitProduction();

// Equation solves and bounds

    /// Solve the omega transport equation
    void solveOmegaEquation(const FaceFluxField& flowRateFace);

    /// Bound omega with viscosity-ratio limit
    void boundOmega();

    /// Solve the k transport equation
    void solveKEquation(const FaceFluxField& flowRateFace);

    /// Bound k to minimum value
    void boundK();

// Post-solve updates and diagnostics

    /// Update SST turbulent viscosity (uses f23_/strainRateMag_)
    void updateTurbulentViscosity();

    /// Update wall shear stress using parallel velocity only
    void updateWallShearStress
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz
    );

    /// Compute normalised k/omega change against the pre-solve snapshots
    void updateResiduals();

    /// Log min/max/mean for k, omega, nut fields
    void logFieldDiagnostics() const;
};
