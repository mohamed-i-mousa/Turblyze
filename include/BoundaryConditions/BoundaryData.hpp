/******************************************************************************
 * @file BoundaryData.h
 * @brief Boundary condition data storage and type management
 * 
 * @class BoundaryData
 * 
 * Stores boundary condition data including type, values, and gradients for
 * individual patches. Supports multiple BC types (FIXED_VALUE, ZERO_GRADIENT,
 * FIXED_GRADIENT, NO_SLIP) with proper scalar/vector data handling.
 * 
 * Key features:
 * - Type-safe storage of boundary condition parameters
 * - Support for scalar and vector boundary values and gradients
 * - Conditional data usage based on boundary condition type
 * - Integration with field-based boundary condition application
 *****************************************************************************/

#ifndef BOUNDARYDATA_H
#define BOUNDARYDATA_H

#include <string>

#include "Scalar.hpp"
#include "Vector.hpp"

/**
 * @enum BCType 
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
 * @enum BCValueType 
 * @brief Enumeration of boundary condition value types
 */
enum class BCValueType 
{
    SCALAR,      ///< Scalar-valued boundary condition
    VECTOR,      ///< Vector-valued boundary condition
    UNDEFINED    ///< Undefined value type
};

class BoundaryData 
{
public:

    /// Default constructor
    BoundaryData() = default;

// Setter methods 

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

// Accesor methods
    
    /** 
     * @brief Get boundary condition type 
     * @return Current BC type 
     */
    BCType type() const { return type_; }
    
    /** 
     * @brief Get value type 
     * @return Type of boundary value (scalar/vector) 
     */
    BCValueType valueType() const { return valueType_; }
    
    /** 
     * @brief Get gradient type 
     * @return Type of boundary gradient (scalar/vector) 
     */
    BCValueType gradientType() const { return gradientType_; }
    
    /** 
     * @brief Get scalar value 
     * @return Current scalar boundary value 
     */
    Scalar scalarValue() const { return scalarValue_; }
    
    /** 
     * @brief Get vector value 
     * @return Current vector boundary value 
     */
    const Vector& vectorValue() const { return vectorValue_; }
    
    /** 
     * @brief Get scalar gradient 
     * @return Current scalar boundary gradient 
     */
    Scalar scalarGradient() const { return scalarGradient_; }
    
    /** 
     * @brief Get vector gradient 
     * @return Current vector boundary gradient 
     */
    const Vector& vectorGradient() const { return vectorGradient_; }

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
    
    /**
     * @brief Get fixed vector gradient
     * @return Fixed vector gradient (normal component)
     * @throws std::runtime_error if not a fixed vector gradient BC
     */
    const Vector& fixedVectorGradient() const;
    
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
    
    /// Vector boundary gradient (normal component)
    Vector vectorGradient_;
};

#endif