/******************************************************************************
 * @file ConvectionScheme.hpp
 * @brief Convection discretization schemes for finite volume method
 *
 * @details This header defines the convection scheme hierarchy for 
 * discretizing convective fluxes in transport equations. The implementation 
 * follows a deferred correction approach where all schemes use stable upwind 
 * matrix coefficients, with higher-order schemes adding explicit correction 
 * terms.
 *
 * @class ConvectionScheme
 *
 * The ConvectionScheme class provides:
 * - Base interface for convective flux discretization
 * - Upwind flux coefficient calculation for matrix assembly
 * - Virtual deferred correction interface for higher-order accuracy
 * - Support for first-order (Upwind), second-order (CDS, SOU) schemes
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <memory>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
#include "CellData.hpp"


class ConvectionScheme
{
public:
    /// Default destructor
    virtual ~ConvectionScheme() = default;

// Public types

    /**
     * @brief Upwind convection flux coefficients for matrix assembly
     *
     * @details
     * All schemes use upwind coefficients in the matrix for stability:
     * - owner = max(F, 0): takes flow when F >= 0
     * - neighbor = min(F, 0): takes flow when F < 0
     */
    struct FluxCoefficients
    {
        Scalar owner;
        Scalar neighbor;
    };

// Public interface

    /**
     * @brief Get upwind convection flux coefficients for matrix assembly
     * @param flowRate Convective volumetric flow rate through the face
     * @return FluxCoefficients with owner and neighbor coefficients
     */
    [[nodiscard("Computed flux coefficients are needed for matrix assembly")]]
    static constexpr FluxCoefficients getFluxCoefficients
    (
        Scalar flowRate
    ) noexcept
    {
        return
        {
            .owner    = std::max(flowRate, S(0.0)),
            .neighbor = std::min(flowRate, S(0.0))
        };
    }

    /**
     * @brief Calculate higher-order deferred correction term
     * @param face Face for interpolation
     * @param phi Cell-centered field values
     * @param gradPhiP Cell gradient at owner cell
     * @param gradPhiN Cell gradient at neighbor cell
     * @param flowRate Volumetric flow rate through face
     * @return Correction flux: flowRate × (φ_highOrder - φ_upwind)
     *
     * Base implementation returns 0 (no correction for first-order upwind).
     * Higher-order schemes override to provide explicit correction terms.
     */
    virtual Scalar calculateCorrection
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const;
};

/**
 * @brief First-order Upwind Differencing Scheme (UDS)
 *
 * @details
 * Discretization: φ_f = φ_upwind (flow-direction dependent)
 *   - F >= 0: φ_f = φ_P (owner cell value)
 *   - F < 0:  φ_f = φ_N (neighbor cell value)
 *
 * Unconditionally stable but diffusive (first-order accurate).
 */
class UpwindScheme final : public ConvectionScheme {};

/**
 * @brief Central Differencing Scheme (CDS) with deferred correction
 *
 * @details
 * Discretization: φ_f = w·φ_P + (1-w)·φ_N (distance-weighted interpolation)
 *
 * Implementation uses deferred correction approach:
 *   - Matrix coefficients: upwind (for stability)
 *   - RHS correction: F × (φ_CDS - φ_upwind)
 *
 * Second-order accurate but can be unstable for high cell Peclet numbers.
 */
class CentralDifferenceScheme final : public ConvectionScheme
{
public:

// Public interface

    /**
     * @brief Calculate deferred correction for CDS
     * @param face Face for interpolation
     * @param phi Cell-centered field values
     * @param gradPhiP Cell gradient at owner (unused for CDS)
     * @param gradPhiN Cell gradient at neighbor (unused for CDS)
     * @param flowRate Volumetric flow rate through face
     * @return Correction flux: F × (φ_CDS - φ_upwind)
     */
    Scalar calculateCorrection
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};

/**
 * @brief Second-order Upwind Scheme (SOU) with gradient reconstruction
 *
 * @details
 * Discretization: φ_f = φ_upwind + ∇φ_upwind · d_{upwind→face}
 *
 * Implementation uses deferred correction approach:
 *   - Matrix coefficients: upwind (for stability)
 *   - RHS correction: F × (φ_SOU - φ_upwind)
 *
 * Balances accuracy (second-order) with stability (upwind bias).
 * More robust than CDS for convection-dominated flows.
 */
class SecondOrderUpwindScheme final : public ConvectionScheme
{
public:

// Public interface

    /**
     * @brief Calculate deferred correction for SOU
     * @param face Face for interpolation
     * @param phi Cell-centered field values
     * @param gradPhiP Cell gradient at owner cell
     * @param gradPhiN Cell gradient at neighbor cell
     * @param flowRate Volumetric flow rate (determines upwind direction)
     * @return Correction flux: F × (φ_SOU - φ_upwind)
     */
    Scalar calculateCorrection
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};

/**
 * @brief Container for equation-specific convection schemes
 *
 * @details
 * Manages convection schemes for different transport equations:
 * - momentum: Ux, Uy, Uz momentum equations
 * - k: Turbulent kinetic energy transport
 * - omega: Specific dissipation rate transport
 * - default: Fallback for any unspecified equation
 *
 * Accessor methods automatically fall back to default scheme when
 * equation-specific scheme is not set.
 */
struct ConvectionSchemes
{
// Scheme storage

    /// Default scheme
    std::unique_ptr<ConvectionScheme> defaultScheme;

    /// Scheme for momentum equations (Ux, Uy, Uz)
    std::unique_ptr<ConvectionScheme> momentumScheme;

    /// Scheme for turbulent kinetic energy (k) equation
    std::unique_ptr<ConvectionScheme> kScheme;

    /// Scheme for specific dissipation rate (omega) equation
    std::unique_ptr<ConvectionScheme> omegaScheme;

// Accessor methods

    /**
     * @brief Get scheme for momentum equations
     * @return Reference to momentum scheme if set, otherwise default scheme
     */
    const ConvectionScheme& momentum() const
    {
        return momentumScheme ? *momentumScheme : *defaultScheme;
    }

    /**
     * @brief Get scheme for k equation
     * @return Reference to k scheme if set, otherwise default scheme
     */
    const ConvectionScheme& k() const
    {
        return kScheme ? *kScheme : *defaultScheme;
    }

    /**
     * @brief Get scheme for omega equation
     * @return Reference to omega scheme if set, otherwise default scheme
     */
    const ConvectionScheme& omega() const
    {
        return omegaScheme ? *omegaScheme : *defaultScheme;
    }
};
