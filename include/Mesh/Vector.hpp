/******************************************************************************
 * @file Vector.h
 * @brief 3D vector class for geometric and mathematical operations in CFD
 * 
 * This header defines a comprehensive 3D vector class that serves as the
 * foundation for all vector-based calculations in the CFD solver. The Vector
 * class provides essential mathematical operations for spatial computations,
 * geometric transformations, and vector algebra required in finite volume
 * discretization and mesh operations.
 * 
 * @class Vector
 * 
 * The Vector class provides:
 * - Component-wise access and manipulation (x, y, z coordinates)
 * - Full arithmetic operations (addition, subtraction, scalar multiplication)
 * - Vector algebra operations (dot product, cross product, normalization)
 * - Geometric calculations (distance, squared distance, magnitude)
 * - Comparison operators for sorting and equality testing
 * - Stream I/O operators for debugging and file output
 * 
 * The class integrates seamlessly with mesh geometry calculations,
 * face normal computations, and field interpolation operations.
 *****************************************************************************/

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include "Scalar.hpp"

class Vector 
{
public:

    /// Default constructor 
    Vector();
    
    /**
     * @brief Constructs vector with specified components
     * @param x_val X component
     * @param y_val Y component 
     * @param z_val Z component
     */
    Vector(Scalar x_val, Scalar y_val, Scalar z_val);

// Setter methods

    /** 
     * @brief Set X component 
     * @param x_val New X component value 
     */
    void setX(Scalar x_val) { x_ = x_val; }
    
    /** 
     * @brief Set Y component 
     * @param y_val New Y component value 
     */
    void setY(Scalar y_val) { y_ = y_val; }
    
    /** 
     * @brief Set Z component 
     * @param z_val New Z component value 
     */
    void setZ(Scalar z_val) { z_ = z_val; }

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

/**
 * @brief Calculates distance between two vectors
 * @param p1 First vector
 * @param p2 Second vector
 * @return Distance between vectors
 */
Scalar distance(const Vector& p1, const Vector& p2);

/**
 * @brief Calculates squared distance between two vectors
 * @param p1 First vector
 * @param p2 Second vector
 * @return Squared distance between vectors
 */
Scalar distanceSquared(const Vector& p1, const Vector& p2);

#endif