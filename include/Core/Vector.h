/**
 * @brief A 3D vector class for geometric operations
 * 
 * Represents a 3D vector with x, y, z components and provides
 * standard vector operations including arithmetic, normalization,
 * and geometric calculations.
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include "Scalar.h"


class Vector 
{
public:
    /// X, Y, Z components of the vector
    Scalar x, y, z;

    /**
     * @brief Default constructor - creates zero vector (0,0,0)
     */
    Vector();
    
    /**
     * @brief Constructs vector with specified components
     * @param x_val X component
     * @param y_val Y component 
     * @param z_val Z component
     */
    Vector(Scalar x_val, Scalar y_val, Scalar z_val);

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
     * @note More efficient than magnitude() for comparisons
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
};

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
 * @note More efficient than distance() for comparisons
 */
Scalar distanceSquared(const Vector& p1, const Vector& p2);

#endif