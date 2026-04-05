/******************************************************************************
 * @file ConvectionScheme.hpp
 * @brief Convection discretization schemes for finite volume method
 *
 * @details This header defines the convection scheme hierarchy for
 * discretizing convective fluxes in transport equations. The implementation
 * follows a deferred correction approach: all schemes use stable first-order
 * upwind coefficients in the implicit matrix, with higher-order schemes
 * adding an explicit correction term to the RHS each iteration.
 *
 * ## Scheme Hierarchy
 *
 * @class ConvectionScheme is the abstract base class. Three concrete schemes
 * are provided:
 *
 * - @class UpwindScheme (1st-order): φ_f = φ_upwind
 *     Selects the upwind cell value based on flow direction.
 *     Unconditionally stable but introduces numerical diffusion.
 *
 * - @class CentralDifferenceScheme (2nd-order): φ_f = w·φ_P + (1-w)·φ_N
 *     Distance-weighted linear interpolation between owner and neighbor.
 *     Second-order accurate but can be unstable at high cell Peclet numbers.
 *
 * - @class SecondOrderUpwindScheme (2nd-order): φ_f = φ_upwind + ∇φ_upwind · d
 *     Gradient-reconstructed face value from the upwind cell center.
 *     Balances accuracy and robustness; preferred for turbulent flows.
 *
 * ## Upwind Matrix Coefficients
 *
 * @struct FluxCoefficients holds the implicit upwind contributions used by
 * all schemes for matrix assembly:
 * - owner    = max(F, 0): active when flow is from owner to neighbor (F >= 0)
 * - neighbor = min(F, 0): active when flow is from neighbor to owner (F < 0)
 *
 * Retrieved via ConvectionScheme::getFluxCoefficients(flowRate).
 *
 * ## ConvectionSchemes Container
 *
 * ConvectionSchemes groups per-equation scheme selections (momentum, k,
 * omega) with a shared default fallback. Accessor methods return the
 * equation-specific scheme when set, otherwise fall back to the default.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <memory>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
#include "CellData.hpp"
#include "ErrorHandler.hpp"


class ConvectionScheme
{
public:
    /// Default destructor
    virtual ~ConvectionScheme() = default;

// Public types

    struct FluxCoefficients
    {
        Scalar owner;
        Scalar neighbor;
    };

// Public interface

    /**
     * @brief Get implicit upwind flux coefficients for matrix assembly
     * @param flowRate Volumetric flow rate through face (F = ρ·U·A)
     * @return FluxCoefficients with owner and neighbor contributions
     *
     * The returned coefficients are used in the implicit matrix assembly:
     * - owner    = max(F, 0): active when flow is from owner to neighbor (F >= 0)
     * - neighbor = min(F, 0): active when flow is from neighbor to owner (F < 0)
     *
     * All schemes use these same upwind coefficients for stability. Higher-order
     * schemes add an explicit correction term to the RHS to achieve improved
     * accuracy while maintaining a stable implicit discretization.
     */
    [[nodiscard("Computed flux coefficients are needed for matrix assembly")]]
    static FluxCoefficients getFluxCoefficients(Scalar flowRate) noexcept
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
     * @return Correction flux: flowRate × (φ_highOrder - φ_upwind).
     *         UpwindScheme returns 0; higher-order schemes return the
     *         explicit correction term added to the RHS.
     */
    
    [[nodiscard("Computed correction is needed for matrix assembly")]]
    virtual Scalar calculateCorrection
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const = 0;
};


class UpwindScheme final : public ConvectionScheme
{
public:

    Scalar calculateCorrection
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};


class CentralDifferenceScheme final : public ConvectionScheme
{
public:

// Public interface

    Scalar calculateCorrection
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};


class SecondOrderUpwindScheme final : public ConvectionScheme
{
public:

// Public interface

    Scalar calculateCorrection
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};

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
        if (!defaultScheme)
        {
            FatalError("Default convection scheme must be set");
        }
        return momentumScheme ? *momentumScheme : *defaultScheme;
    }

    /**
     * @brief Get scheme for k equation
     * @return Reference to k scheme if set, otherwise default scheme
     */
    const ConvectionScheme& k() const
    {
        if (!defaultScheme)
        {
            FatalError("Default convection scheme must be set");
        }
        return kScheme ? *kScheme : *defaultScheme;
    }

    /**
     * @brief Get scheme for omega equation
     * @return Reference to omega scheme if set, otherwise default scheme
     */
    const ConvectionScheme& omega() const
    {
        if (!defaultScheme)
        {
            FatalError("Default convection scheme must be set");
        }
        return omegaScheme ? *omegaScheme : *defaultScheme;
    }
};
