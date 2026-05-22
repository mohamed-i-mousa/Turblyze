/******************************************************************************
 * @file ConvectionScheme.h
 * @brief Per-equation convection scheme container
 *
 * @class ConvectionScheme groups per-equation scheme selections (momentum, k,
 * omega) with a shared default fallback. Accessor methods return the
 * equation-specific scheme when set, otherwise fall back to the default.
 ******************************************************************************/

#pragma once

#include <memory>

#include "ConvectionSchemes.h"


class ConvectionScheme
{
public:

    ConvectionScheme() = default;

    /// Copy constructor and assignment - Not copyable (unique_ptr members)
    ConvectionScheme(const ConvectionScheme&) = delete;
    ConvectionScheme& operator=(const ConvectionScheme&) = delete;

    /// Move constructor and assignment
    ConvectionScheme(ConvectionScheme&&) = default;
    ConvectionScheme& operator=(ConvectionScheme&&) = default;

    /// Destructor
    ~ConvectionScheme() noexcept = default;

// Scheme storage

    /// Default scheme
    std::unique_ptr<ConvectionSchemes> defaultScheme;

    /// Scheme for momentum equations (Ux, Uy, Uz)
    std::unique_ptr<ConvectionSchemes> momentumScheme;

    /// Scheme for turbulent kinetic energy (k) equation
    std::unique_ptr<ConvectionSchemes> kScheme;

    /// Scheme for specific dissipation rate (omega) equation
    std::unique_ptr<ConvectionSchemes> omegaScheme;

// Accessor methods

    /// Get scheme for momentum equations
    const ConvectionSchemes& momentum() const
    {
        return resolve(momentumScheme.get());
    }

    /// Get scheme for k equation
    const ConvectionSchemes& k() const
    {
        return resolve(kScheme.get());
    }

    /// Get scheme for omega equation
    const ConvectionSchemes& omega() const
    {
        return resolve(omegaScheme.get());
    }

private:

    const ConvectionSchemes& resolve
    (
        const ConvectionSchemes* specific
    ) const
    {
        if (specific) return *specific;
        if (!defaultScheme)
        {
            FatalError("Default convection scheme must be set");
        }
        return *defaultScheme;
    }
};
