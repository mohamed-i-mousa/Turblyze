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
 * - Enumerations for supported BC types (FIXED_VALUE, ZERO_GRADIENT, etc.)
 * - Validation logic to ensure data consistency
 * - Conditional accessors based on the active configuration
 *****************************************************************************/

#pragma once

#include "Scalar.h"


enum class BCType
{
    FIXED_VALUE,            ///< Fixed value (Dirichlet) boundary condition
    FIXED_GRADIENT,         ///< Fixed gradient (Neumann) boundary condition
    ZERO_GRADIENT,          ///< Zero gradient boundary condition
    NO_SLIP,                ///< No-slip wall boundary condition
    K_WALL_FUNCTION,        ///< kWallFunction-like boundary condition
    OMEGA_WALL_FUNCTION,    ///< omegaWallFunction-like boundary condition
    NUT_WALL_FUNCTION,      ///< nutWallFunction-like boundary condition
    UNDEFINED               ///< Undefined boundary condition type
};


class BoundaryData
{
public:

// Setter methods

    /**
     * @brief Set fixed scalar value boundary condition
     * @param scalarValue Scalar value to fix at boundary
     */
    void setFixedValue(Scalar scalarValue) noexcept;

    /**
     * @brief Set fixed scalar gradient boundary condition
     * @param scalarGradient Scalar gradient (normal component) at boundary
     */
    void setFixedGradient(Scalar scalarGradient) noexcept;

    /// Set zero gradient boundary condition
    void setZeroGradient() noexcept;

    /// Set no-slip boundary condition (for velocity)
    void setNoSlip() noexcept;

    /**
     * @brief Set wall function boundary condition type
     * @param wallType The wall function type (K_WALL_FUNCTION, OMEGA_WALL_FUNCTION, or NUT_WALL_FUNCTION)
     */
    void setWallFunctionType(BCType wallType) noexcept;

// Accessor methods

    /**
     * @brief Get boundary condition type
     * @return Current BC type
     */
    [[nodiscard]] BCType type() const noexcept { return type_; }

    /**
     * @brief Get scalar value
     * @return Current scalar boundary value
     */
    [[nodiscard]] Scalar scalarValue() const noexcept { return scalarValue_; }

    /**
     * @brief Get scalar gradient
     * @return Current scalar boundary gradient
     */
    [[nodiscard]] Scalar scalarGradient() const noexcept
    {
        return scalarGradient_;
    }

    /**
     * @brief Get fixed scalar value
     * @return Fixed scalar value
     * @note Terminates the program if not a fixed scalar value BC
     */
    [[nodiscard]] Scalar fixedScalarValue() const;

    /**
     * @brief Get fixed scalar gradient
     * @return Fixed scalar gradient (normal component)
     * @note Terminates the program if not a fixed scalar gradient BC
     */
    [[nodiscard]] Scalar fixedScalarGradient() const;

private:

// Private members

    /// Boundary condition type
    BCType type_ = BCType::UNDEFINED;

    /// Scalar boundary value
    Scalar scalarValue_ = S(0.0);

    /// Scalar boundary gradient (normal component)
    Scalar scalarGradient_ = S(0.0);
};
