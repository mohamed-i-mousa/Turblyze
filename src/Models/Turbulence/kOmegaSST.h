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

// Standard library headers
#include <cmath>
#include <vector>
#include <limits>

// Project headers
#include "RANS.h"
#include "Scalar.h"
#include "Vector.h"
#include "Mesh.h"
#include "CellData.h"
#include "FaceData.h"

// ****************************** class kOmegaSST *****************************

class kOmegaSST final : public RANS
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
        Scalar nu,
        Scalar initialK,
        Scalar initialOmega,
        Scalar alphaK,
        Scalar alphaOmega,
        bool debug
    );

    /// Copy constructor and assignment - Not copyable (non-copyable base)
    kOmegaSST(const kOmegaSST&) = delete;
    kOmegaSST& operator=(const kOmegaSST&) = delete;

    /// Move constructor and assignment - Not movable (non-copyable base)
    kOmegaSST(kOmegaSST&&) = delete;
    kOmegaSST& operator=(kOmegaSST&&) = delete;

    /// Destructor
    ~kOmegaSST() noexcept override;

// ****************************** Solve kOmegaSST *****************************

    /// Solve turbulence equations for current iteration
    void solve
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz,
        const FaceFluxField& flowRateFace,
        const TensorField& gradU
    ) override;

// ***************************** Accessor Methods *****************************

    /// Get specific dissipation rate field (native k-omega accessor)
    [[nodiscard]] const ScalarField& omega() const noexcept
    {
        return omega_;
    }

    /// Dissipation field for the base contract (omega for k-omega SST)
    [[nodiscard]] const ScalarField& dissipation() const noexcept override
    {
        return omega();
    }

    /// Dissipation field name for output labelling
    [[nodiscard]] NameRef dissipationName() const noexcept override
    {
        return "omega";
    }

    /// This model carries turbulence
    [[nodiscard]] bool isTurbulent() const noexcept override
    {
        return true;
    }

// ****************************** Model Constants *****************************

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

// ************************ Inlet Condition Calculators ***********************

    /// Calculate inlet/initial turbulent kinetic energy
    [[nodiscard]] static Scalar inletK
    (
        const Vector& velocity,
        Scalar turbulenceIntensity
    ) noexcept
    {
        return RANS::inletK(velocity, turbulenceIntensity);
    }

    /// Calculate inlet/initial specific dissipation rate
    [[nodiscard]] static Scalar inletOmega
    (
        Scalar k,
        Scalar hydraulicDiameter
    ) noexcept;

// ****************************** Private Members *****************************

private:

// SST-specific parameters

    /// Optional SST F3 switch
    bool useF3_ = false;

    /// Maximum turbulent-to-laminar viscosity ratio for omega bound
    Scalar maxViscosityRatio_ = S(1e5);

// Turbulence fields

    /// Specific dissipation rate
    ScalarField omega_{S(1.0)};

// Residual tracking

    /// Previous-iteration omega snapshot for residual computation
    ScalarField omegaPrev_;

    /// Face diffusivity for the k equation
    FaceFluxField gammaKFace_;

    /// Face diffusivity for the omega equation
    FaceFluxField gammaOmegaFace_;

// Wall distance and wall-function state

    /// Dynamic omega wall-function values on faces
    FaceData<Scalar> omegaWall_
    {
        std::numeric_limits<Scalar>::quiet_NaN()
    };

    /// Area-weighted wall-function omega per wall cell
    ScalarList wallCellOmega_;

    /// Precomputed Cmu^0.25 (avoids repeated std::pow calls)
    const Scalar Cmu25_ = std::sqrt(std::sqrt(coeffs_.betaStar));

// ****************************** Private Methods *****************************

private:

// Utility helpers

    /// SST wall-function Cμ^0.25, taken as β*^0.25
    [[nodiscard]] Scalar cmu25() const noexcept override
    {
        return Cmu25_;
    }

    /// Blend two constants using SST blending function
    [[nodiscard]] static Scalar blend
    (
        Scalar f,
        Scalar cInner,
        Scalar cOuter
    ) noexcept
    {
        return f * (cInner - cOuter) + cOuter;
    }

// Initialization helpers

    /// Initialize turbulent viscosity with k/omega estimate
    void initializeTurbulentViscosity();

// Wall-function helpers

    /// Update wall-function nut on wall faces (nutkWallFunction)
    void updateNutWall();

    /// Update dynamic omega wall-function boundary values per face
    void updateOmegaWallValues();

    /// Pre-set wall-cell omega to area-weighted wall-function values
    void applyOmegaWallCellValues();

    /// Override k production at wall-adjacent cells
    void overrideWallCellProduction
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz,
        ScalarField& Pk
    );

// Per-iteration model terms

    /// Compute k production term
    [[nodiscard]] ScalarField kProduction
    (
        const ScalarField& strainRateMag
    ) const;

    /// Compute cross-diffusion term
    [[nodiscard]] ScalarField crossDiffusion
    (
        const VectorField& gradOmega
    ) const;

    /// Compute F1 blending function
    [[nodiscard]] ScalarField blendingF1
    (
        const ScalarField& CDkOmega
    ) const;

    /// Compute F2 blending function
    [[nodiscard]] ScalarField blendingF2() const;

    /// Compute F3 blending switch function
    [[nodiscard]] ScalarField blendingF3() const;

    /// Compute F23 blending selector
    [[nodiscard]] ScalarField blendingF23
    (
        const ScalarField& f2,
        const ScalarField& f3
    ) const;

    /// Compute omega production
    [[nodiscard]] ScalarField omegaProduction
    (
        const ScalarField& f1,
        const ScalarField& strainRateMag
    ) const;

    /// Compute k equation effective diffusivity
    [[nodiscard]] ScalarField computeGammaK(const ScalarField& f1) const;

    /// Compute omega equation effective diffusivity
    [[nodiscard]] ScalarField computeGammaOmega(const ScalarField& f1) const;

    /// Limit Pk/POmega in place using f1/f23/strainRateMag
    void limitProduction
    (
        const ScalarField& f1,
        const ScalarField& f23,
        const ScalarField& strainRateMag,
        ScalarField& Pk,
        ScalarField& POmega
    ) const;

// Equation solves and bounds

    /// Solve the omega transport equation
    void solveOmegaEquation
    (
        const FaceFluxField& flowRateFace,
        const ScalarField& divU,
        const ScalarField& f1,
        const ScalarField& CDkOmega,
        const ScalarField& POmega,
        const VectorField& gradOmega
    );

    /// Bound omega with viscosity-ratio limit
    void boundOmega();

    /// Solve the k transport equation
    void solveKEquation
    (
        const FaceFluxField& flowRateFace,
        const ScalarField& divU,
        const ScalarField& f1,
        const ScalarField& Pk
    );

    /// Bound k to minimum value
    void boundK();

// Post-solve updates and diagnostics

    /// Compute SST turbulent viscosity
    [[nodiscard]] ScalarField computeTurbulentViscosity
    (
        const ScalarField& f23,
        const ScalarField& strainRateMag
    ) const;

    /// Log min/max/mean for k, omega, nut fields
    void logFieldDiagnostics() const;
};
