/******************************************************************************
 * @file ConvectionSchemes.h
 * @brief Abstract base class and concrete schemes for convection discretization
 *
 * @details This header defines the convection scheme hierarchy for
 * discretizing convective fluxes in transport equations. The implementation
 * follows a deferred correction approach: all schemes use stable first-order
 * upwind coefficients in the implicit matrix, with higher-order schemes
 * adding an explicit correction term to the RHS each iteration.
 *
 *
 * @class ConvectionSchemes is the abstract base class. Three concrete schemes
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
 *
 * @struct FluxCoefficients holds the implicit upwind contributions used by
 * all schemes for matrix assembly:
 * - owner    = max(F, 0): active when flow is from owner to neighbor (F >= 0)
 * - neighbor = min(F, 0): active when flow is from neighbor to owner (F < 0)
 *
 *****************************************************************************/

#pragma once

#include <algorithm>

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "CellData.h"
#include "ErrorHandler.h"


class ConvectionSchemes
{
public:

    /// Copy constructor and assignment - Slicing problem
    ConvectionSchemes(const ConvectionSchemes&) = delete;
    ConvectionSchemes& operator=(const ConvectionSchemes&) = delete;

    /// Move constructor and assignment - Slicing problem
    ConvectionSchemes(ConvectionSchemes&&) = delete;
    ConvectionSchemes& operator=(ConvectionSchemes&&) = delete;

    /// Default destructor
    virtual ~ConvectionSchemes() = default;

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
    [[nodiscard]] static FluxCoefficients getFluxCoefficients
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
     * @return Correction flux: flowRate × (φ_highOrder - φ_upwind).
     *         UpwindScheme returns 0; higher-order schemes return the
     *         explicit correction term added to the RHS.
     */
    [[nodiscard]] virtual Scalar correction
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const = 0;

protected:

    /// Default constructor
    ConvectionSchemes() = default;
};


class UpwindScheme final : public ConvectionSchemes
{
public:

// Public interface

    Scalar correction
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};


class CentralDifferenceScheme final : public ConvectionSchemes
{
public:

// Public interface

    Scalar correction
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};


class SecondOrderUpwindScheme final : public ConvectionSchemes
{
public:

// Public interface

    Scalar correction
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};
