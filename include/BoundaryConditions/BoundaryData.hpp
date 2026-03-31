/******************************************************************************
 * @file BoundaryData.hpp
 * @brief Boundary condition data storage and type management
 *
 * @details This header defines the BoundaryData class, which acts as a data
 * transfer object holding the specific parameters for boundary conditions
 * (e.g., "velocity on inlet patch").
 *
 * @enum BCType
 * - Enumeration of boundary condition types
 * 
 * @enum BCValueType
 * - Enumeration of boundary condition value types
 * 
 * @class BoundaryData
 * - Type-safe storage for BC parameters (scalar/vector values, gradients)
 * - Enumerations for supported BC types (FIXED_VALUE, ZERO_GRADIENT, etc.)
 * - Validation logic to ensure data consistency
 * - Conditional accessors based on the active configuration
 *****************************************************************************/

#pragma once

#include "Scalar.hpp"
#include "Vector.hpp"


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


enum class BCValueType
{
    SCALAR,                 ///< Scalar-valued boundary condition
    VECTOR,                 ///< Vector-valued boundary condition
    UNDEFINED               ///< Undefined value type
};


class BoundaryData
{
public:

    /// Default constructor
    BoundaryData() noexcept = default;

// Setter methods

    /**
     * @brief Set fixed scalar value boundary condition
     * @param scalarValue Scalar value to fix at boundary
     */
    void setFixedValue(Scalar scalarValue) noexcept;

    /**
     * @brief Set fixed vector value boundary condition
     * @param vectorValue Vector value to fix at boundary
     */
    void setFixedValue(const Vector& vectorValue) noexcept;

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
    BCType type() const noexcept { return type_; }

    /**
     * @brief Get value type
     * @return Type of boundary value (scalar/vector)
     */
    BCValueType valueType() const noexcept { return valueType_; }

    /**
     * @brief Get gradient type
     * @return Type of boundary gradient (scalar/vector)
     */
    BCValueType gradientType() const noexcept { return gradientType_; }

    /**
     * @brief Get scalar value
     * @return Current scalar boundary value
     */
    Scalar scalarValue() const noexcept { return scalarValue_; }

    /**
     * @brief Get vector value
     * @return Current vector boundary value
     */
    const Vector& vectorValue() const noexcept { return vectorValue_; }

    /**
     * @brief Get scalar gradient
     * @return Current scalar boundary gradient
     */
    Scalar scalarGradient() const noexcept { return scalarGradient_; }

    /**
     * @brief Get fixed scalar value
     * @return Fixed scalar value
     * @throws std::runtime_error if not a fixed scalar value BC
     */
    Scalar fixedScalarValue() const;

    /**
     * @brief Get fixed vector value
     * @return Fixed vector value
     * @throws std::runtime_error if not a fixed vector value BC
     */
    const Vector& fixedVectorValue() const;

    /**
     * @brief Get fixed scalar gradient
     * @return Fixed scalar gradient (normal component)
     * @throws std::runtime_error if not a fixed scalar gradient BC
     */
    Scalar fixedScalarGradient() const;

private:

// Private members

    /// Boundary condition type
    BCType type_ = BCType::UNDEFINED;

    /// Type of boundary value (scalar or vector)
    BCValueType valueType_ = BCValueType::UNDEFINED;

    /// Scalar boundary value
    Scalar scalarValue_ = S(0.0);

    /// Vector boundary value
    Vector vectorValue_;

    /// Type of boundary gradient (scalar or vector)
    BCValueType gradientType_ = BCValueType::UNDEFINED;

    /// Scalar boundary gradient (normal component)
    Scalar scalarGradient_ = S(0.0);
};
