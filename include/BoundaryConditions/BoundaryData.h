#ifndef BOUNDARYDATA_H
#define BOUNDARYDATA_H

#include <string>

#include "Scalar.h"
#include "Vector.h"

/**
 * @brief Enumeration of boundary condition types
 */
enum class BCType 
{
    FIXED_VALUE,     ///< Fixed value (Dirichlet) boundary condition
    FIXED_GRADIENT,  ///< Fixed gradient (Neumann) boundary condition
    ZERO_GRADIENT,   ///< Zero gradient boundary condition
    NO_SLIP,         ///< No-slip wall boundary condition
    UNDEFINED        ///< Undefined boundary condition type
};

/**
 * @brief Enumeration of boundary condition value types
 */
enum class BCValueType 
{
    SCALAR,      ///< Scalar-valued boundary condition
    VECTOR,      ///< Vector-valued boundary condition
    UNDEFINED    ///< Undefined value type
};

/**
 * @brief Stores boundary condition data for a patch
 * 
 * This struct encapsulates boundary condition information including type,
 * values, and gradients. Depending on the boundary condition type, either
 * scalar or vector values/gradients are used while others remain undefined.
 * 
 * For gradient boundary conditions, only the normal component of the
 * gradient ∂φ/∂n is specified.
 */
struct BoundaryData 
{
    /// Boundary condition type
    BCType type = BCType::UNDEFINED;

    /// Type of boundary value (scalar or vector)
    BCValueType valueType = BCValueType::UNDEFINED;
    
    /// Scalar boundary value
    Scalar scalarValue = S(0.0);
    
    /// Vector boundary value
    Vector vectorValue; 

    /// Type of boundary gradient (scalar or vector)
    BCValueType gradientType = BCValueType::UNDEFINED;
    
    /// Scalar boundary gradient (normal component)
    Scalar scalarGradient = S(0.0);
    
    /// Vector boundary gradient (normal component)
    Vector vectorGradient;         

    /**
     * @brief Default constructor
     */
    BoundaryData() = default;

    /**
     * @brief Set fixed scalar value boundary condition
     * @param s_val Scalar value to fix at boundary
     */
    void setFixedValue(Scalar s_val);
    
    /**
     * @brief Set fixed vector value boundary condition
     * @param v_val Vector value to fix at boundary
     */
    void setFixedValue(const Vector& v_val);
    
    /**
     * @brief Set fixed scalar gradient boundary condition
     * @param s_grad Scalar gradient (normal component) at boundary
     */
    void setFixedGradient(Scalar s_grad);
    
    /**
     * @brief Set fixed vector gradient boundary condition
     * @param v_grad Vector gradient (normal component) at boundary
     */
    void setFixedGradient(const Vector& v_grad);
    
    /**
     * @brief Set zero gradient boundary condition
     */
    void setZeroGradient();
    
    
    /**
     * @brief Set no-slip boundary condition (for velocity)
     */
    void setNoSlip();

    /**
     * @brief Get fixed scalar value
     * @return Fixed scalar value
     * @throws std::runtime_error if not a fixed scalar value BC
     */
    Scalar getFixedScalarValue() const;
    
    /**
     * @brief Get fixed vector value
     * @return Fixed vector value
     * @throws std::runtime_error if not a fixed vector value BC
     */
    const Vector& getFixedVectorValue() const;
    
    /**
     * @brief Get fixed scalar gradient
     * @return Fixed scalar gradient (normal component)
     * @throws std::runtime_error if not a fixed scalar gradient BC
     */
    Scalar getFixedScalarGradient() const;
    
    /**
     * @brief Get fixed vector gradient
     * @return Fixed vector gradient (normal component)
     * @throws std::runtime_error if not a fixed vector gradient BC
     */
    const Vector& getFixedVectorGradient() const;
};

#endif