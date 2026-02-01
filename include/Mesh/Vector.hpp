/******************************************************************************
 * @file Vector.hpp
 * @brief 3D vector class for geometric and mathematical operations in CFD
 * 
 * This header defines a 3D vector class that serves as the foundation for 
 * all vector-based calculations in the CFD solver. The Vector class provides
 * essential mathematical operations required in the finite volume
 * discretization and mesh operations.
 * 
 * @class Vector
 * 
 * The Vector class provides:
 * - Components access and manipulation (x, y, z coordinates)
 * - Arithmetic operations (addition, subtraction, scalar multiplication)
 * - Vector operations (dot product, cross product, normalization)
 * - Comparison operators for sorting and equality testing
 * - Stream I/O operators for debugging
 *****************************************************************************/

#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <iostream>
#include "Scalar.hpp"


class Vector 
{
public:

    /// Default constructor 
    Vector();
    
    /**
     * @brief Constructs vector with specified components
     * @param x_value X component
     * @param y_value Y component 
     * @param z_value Z component
     */
    Vector(Scalar x_value, Scalar y_value, Scalar z_value);

// Setter methods

    /** 
     * @brief Set X component 
     * @param x_value New X component value 
     */
    void setX(Scalar x_value) { x_ = x_value; }
    
    /** 
     * @brief Set Y component 
     * @param y_value New Y component value 
     */
    void setY(Scalar y_value) { y_ = y_value; }
    
    /** 
     * @brief Set Z component 
     * @param z_value New Z component value 
     */
    void setZ(Scalar z_value) { z_ = z_value; }

// Accessor methods

    /** 
     * @brief Get X component 
     * @return X component value 
     */
    Scalar x() const { return x_; }
    
    /** 
     * @brief Get Y component 
     * @return Y component value 
     */
    Scalar y() const { return y_; }
    
    /** 
     * @brief Get Z component 
     * @return Z component value 
     */
    Scalar z() const { return z_; }

// Operator methods

    /**
     * @brief Vector addition operator
     * @param other Vector to add
     * @return Sum of the two vectors
     */
    Vector operator+(const Vector& other) const;
    
    /**
     * @brief Vector subtraction operator
     * @param other Vector to subtract
     * @return Difference of the two vectors
     */
    Vector operator-(const Vector& other) const;
    
    /**
     * @brief Scalar multiplication operator
     * @param scalar Scalar to multiply by
     * @return Vector scaled by scalar
     */
    Vector operator*(Scalar scalar) const;
    
    /**
     * @brief Scalar division operator
     * @param scalar Scalar to divide by
     * @return Vector divided by scalar
     * @throws std::runtime_error if scalar is near zero
     */
    Vector operator/(Scalar scalar) const;
    
    /**
     * @brief Compound addition assignment operator
     * @param other Vector to add
     * @return Reference to this vector
     */
    Vector& operator+=(const Vector& other);
    
    /**
     * @brief Compound subtraction assignment operator
     * @param other Vector to subtract
     * @return Reference to this vector
     */
    Vector& operator-=(const Vector& other);
    
    /**
     * @brief Compound multiplication assignment operator
     * @param scalar Scalar to multiply by
     * @return Reference to this vector
     */
    Vector& operator*=(Scalar scalar);
    
    /**
     * @brief Compound division assignment operator
     * @param scalar Scalar to divide by
     * @return Reference to this vector
     * @throws std::runtime_error if scalar is near zero
     */
    Vector& operator/=(Scalar scalar);
    
    /**
     * @brief Equality comparison operator
     * @param other Vector to compare with
     * @return True if vectors are equal within tolerance
     */
    bool operator==(const Vector& other) const;
    
    /**
     * @brief Inequality comparison operator
     * @param other Vector to compare with
     * @return True if vectors are not equal
     */
    bool operator!=(const Vector& other) const;

// Geometric calculations

    /**
     * @brief Calculates squared magnitude of vector
     * @return Squared magnitude (x² + y² + z²)
     */
    Scalar magnitudeSquared() const;
    
    /**
     * @brief Calculates magnitude (length) of vector
     * @return Vector magnitude
     */
    Scalar magnitude() const;
    
    /**
     * @brief Normalizes this vector to unit length
     * @return Reference to this vector for chaining
     * @throws std::runtime_error if vector has zero magnitude
     */
    Vector& normalize();
    
    /**
     * @brief Returns normalized copy of this vector
     * @return Normalized vector
     * @throws std::runtime_error if vector has zero magnitude
     */
    Vector normalized() const;

    /// friend function for operator<<
    friend std::ostream& operator<<(std::ostream& os, const Vector& p);

private:

    /// X, Y, Z components of the vector 
    Scalar x_, y_, z_;
};

// Non-member methods 

/**
 * @brief Scalar multiplication operator (scalar * vector)
 * @param scalar Scalar multiplier
 * @param p Vector to multiply
 * @return Scaled vector
 */
Vector operator*(Scalar scalar, const Vector& p);

/**
 * @brief Stream output operator for Vector
 * @param os Output stream
 * @param p Vector to output
 * @return Reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Vector& p);

/**
 * @brief Computes dot product of two vectors
 * @param p1 First vector
 * @param p2 Second vector
 * @return Scalar dot product
 */
Scalar dot(const Vector& p1, const Vector& p2);

/**
 * @brief Computes cross product of two vectors
 * @param p1 First vector
 * @param p2 Second vector
 * @return Cross product vector
 */
Vector cross(const Vector& p1, const Vector& p2);

#endif // VECTOR_HPP