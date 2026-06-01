/******************************************************************************
 * @file ConvectionSchemes.h
 * @brief Abstract base class for convection discretization schemes
 *
 * @details This header defines the polymorphic convection scheme interface for
 * discretizing convective fluxes in transport equations. All schemes use
 * stable first-order upwind coefficients in the implicit matrix. Higher-order
 * schemes add an explicit deferred correction term through correction().
 *
 * @struct FluxCoefficients holds the implicit upwind contributions used by
 * all schemes for matrix assembly:
 * - owner    = max(F, 0): active when flow is from owner to neighbor (F >= 0)
 * - neighbor = min(F, 0): active when flow is from neighbor to owner (F < 0)
 *
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <algorithm>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "CellData.h"

// ************************** class ConvectionSchemes *************************

class ConvectionSchemes
{
public:

// ************************* Special Member Functions *************************

    /// Copy constructor and assignment - Slicing problem
    ConvectionSchemes(const ConvectionSchemes&) = delete;
    ConvectionSchemes& operator=(const ConvectionSchemes&) = delete;

    /// Move constructor and assignment - Slicing problem
    ConvectionSchemes(ConvectionSchemes&&) = delete;
    ConvectionSchemes& operator=(ConvectionSchemes&&) = delete;

    /// Default destructor
    virtual ~ConvectionSchemes() = default;

// ************************** struct FluxCoefficients *************************

    struct FluxCoefficients
    {
        Scalar owner;
        Scalar neighbor;
    };

// ****************************** Public Methods ******************************

    /// Get implicit upwind flux coefficients for matrix assembly
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

    /// Calculate higher-order deferred correction term
    [[nodiscard]] virtual Scalar correction
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const = 0;

// ***************************** Protected Methods ****************************

protected:

    /// Default constructor
    ConvectionSchemes() = default;
};
