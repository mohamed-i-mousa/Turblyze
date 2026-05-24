/******************************************************************************
 * @file BoundaryData.h
 * @brief Boundary condition data storage and type management
 *
 * @details This header defines the BoundaryData class, which acts as a data
 * transfer object holding the specific parameters for boundary conditions
 * (e.g., "velocity on inlet patch").
 *
 * @enum BCType
 * - Enumeration of boundary condition types
 *
 * @class BoundaryData
 * - Type-safe storage for BC parameters (scalar values, gradients)
 * - Enumerations for supported BC types (fixedValue, zeroGradient, etc.)
 * - Validation logic to ensure data consistency
 * - Conditional accessors based on the active configuration
 *****************************************************************************/

#pragma once

#include <string_view>

#include "Scalar.h"


enum class BCType
{
    fixedValue,             ///< Fixed value (Dirichlet) boundary condition
    fixedGradient,          ///< Fixed gradient (Neumann) boundary condition
    zeroGradient,           ///< Zero gradient boundary condition
    noSlip,                 ///< No-slip wall boundary condition
    kWallFunction,          ///< kWallFunction-like boundary condition
    omegaWallFunction,      ///< omegaWallFunction-like boundary condition
    nutWallFunction,        ///< nutWallFunction-like boundary condition
    undefined               ///< Undefined boundary condition type
};


class BoundaryData
{
public:

// Setter methods

    /// Set fixed scalar value boundary condition
    void setFixedValue(Scalar scalarValue) noexcept;

    /// Set fixed scalar gradient boundary condition
    void setFixedGradient(Scalar scalarGradient) noexcept;

    /// Set zero gradient boundary condition
    void setZeroGradient() noexcept;

    /// Set no-slip boundary condition (for velocity)
    void setNoSlip() noexcept;

    /// Set wall function boundary condition type
    void setWallFunctionType(BCType wallType) noexcept;

// Accessor methods

    /// Get boundary condition type
    [[nodiscard]] BCType type() const noexcept { return type_; }

    /// Get fixed scalar value (fixedValue or noSlip)
    [[nodiscard]] Scalar fixedScalarValue() const noexcept;

    /// Get fixed scalar gradient (fixedGradient)
    [[nodiscard]] Scalar fixedScalarGradient() const noexcept;

private:

// Private members

    /// Boundary condition type
    BCType type_ = BCType::undefined;

    /// Scalar boundary value
    Scalar scalarValue_ = S(0.0);

    /// Scalar boundary gradient (normal component)
    Scalar scalarGradient_ = S(0.0);
};

// Non-member methods

/// Convert BCType to human-readable string
[[nodiscard]] std::string_view bcTypeToString(BCType bctype) noexcept;
