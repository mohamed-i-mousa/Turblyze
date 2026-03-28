/******************************************************************************
 * @file Vector.hpp
 * @brief 3D vector class for geometric and mathematical operations in CFD
 *
 * @details This header defines a 3D vector class that serves as the foundation
 * for all vector-based calculations in the CFD solver. The Vector class
 * provides essential mathematical operations required in the finite volume
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

#pragma once

#include <cmath>
#include <iosfwd>
#include <stdexcept>

#include "Scalar.hpp"

class Vector
{
public:

    /// Default constructor
    Vector() noexcept = default;

    /**
     * @brief Constructs vector with specified components
     * @param xValue X component
     * @param yValue Y component
     * @param zValue Z component
     */
    Vector(Scalar xValue, Scalar yValue, Scalar zValue) noexcept
        : x_(xValue), y_(yValue), z_(zValue) {}

// Setter methods

    /**
     * @brief Set X component
     * @param xValue New X component value
     */
    void setX(Scalar xValue) noexcept { x_ = xValue; }

    /**
     * @brief Set Y component
     * @param yValue New Y component value
     */
    void setY(Scalar yValue) noexcept { y_ = yValue; }

    /**
     * @brief Set Z component
     * @param zValue New Z component value
     */
    void setZ(Scalar zValue) noexcept { z_ = zValue; }

// Accessor methods

    /**
     * @brief Get X component
     * @return X component value
     */
    Scalar x() const noexcept { return x_; }

    /**
     * @brief Get Y component
     * @return Y component value
     */
    Scalar y() const noexcept { return y_; }

    /**
     * @brief Get Z component
     * @return Z component value
     */
    Scalar z() const noexcept { return z_; }

// Operator methods

    /**
     * @brief Vector addition operator
     * @param other Vector to add
     * @return Sum of the two vectors
     */
    Vector operator+(const Vector& other) const noexcept
    {
        Vector result(*this);
        result += other;
        return result;
    }

    /**
     * @brief Vector subtraction operator
     * @param other Vector to subtract
     * @return Difference of the two vectors
     */
    Vector operator-(const Vector& other) const noexcept
    {
        Vector result(*this);
        result -= other;
        return result;
    }

    /**
     * @brief Scalar multiplication operator
     * @param scalar Scalar to multiply by
     * @return Vector scaled by scalar
     */
    Vector operator*(Scalar scalar) const noexcept
    {
        Vector result(*this);
        result *= scalar;
        return result;
    }

    /**
     * @brief Scalar division operator
     * @param scalar Scalar to divide by
     * @return Vector divided by scalar
     * @throws std::runtime_error if scalar is near zero
     */
    Vector operator/(Scalar scalar) const
    {
        Vector result(*this);
        result /= scalar;
        return result;
    }

    /**
     * @brief Compound addition assignment operator
     * @param other Vector to add
     * @return Reference to this vector
     */
    Vector& operator+=(const Vector& other) noexcept
    {
        x_ += other.x_;
        y_ += other.y_;
        z_ += other.z_;

        return *this;
    }

    /**
     * @brief Compound subtraction assignment operator
     * @param other Vector to subtract
     * @return Reference to this vector
     */
    Vector& operator-=(const Vector& other) noexcept
    {
        x_ -= other.x_;
        y_ -= other.y_;
        z_ -= other.z_;

        return *this;
    }

    /**
     * @brief Compound multiplication assignment operator
     * @param scalar Scalar to multiply by
     * @return Reference to this vector
     */
    Vector& operator*=(Scalar scalar) noexcept
    {
        x_ *= scalar;
        y_ *= scalar;
        z_ *= scalar;

        return *this;
    }

    /**
     * @brief Compound division assignment operator
     * @param scalar Scalar to divide by
     * @return Reference to this vector
     * @throws std::runtime_error if scalar is near zero
     */
    Vector& operator/=(Scalar scalar)
    {
        if (std::abs(scalar) <= vSmallValue)
        {
            throw
                std::runtime_error
                (
                    "Error: Division by zero in Vector::operator/="
                );
        }

        Scalar inv = S(1.0) / scalar;
        x_ *= inv;
        y_ *= inv;
        z_ *= inv;

        return *this;
    }

    /**
     * @brief Equality comparison operator
     * @param other Vector to compare with
     * @return True if vectors are equal within tolerance
     * @note the compiler will auto-synthesize operator!= from operator==
     */
    bool operator==(const Vector& other) const noexcept
    {
        return (std::abs(x_ - other.x_) <= vSmallValue)
            && (std::abs(y_ - other.y_) <= vSmallValue)
            && (std::abs(z_ - other.z_) <= vSmallValue);
    }

// Geometric calculations

    /**
     * @brief Calculates squared magnitude of vector
     * @return Squared magnitude (x² + y² + z²)
     */
    Scalar magnitudeSquared() const noexcept
    {
        return x_ * x_ + y_ * y_ + z_ * z_;
    }

    /**
     * @brief Calculates magnitude (length) of vector
     * @return Vector magnitude
     */
    Scalar magnitude() const noexcept
    {
        return std::sqrt(magnitudeSquared());
    }

    /**
     * @brief Returns normalized copy of this vector
     * @return Normalized vector
     * @throws std::runtime_error if vector has zero magnitude
     */
    Vector normalized() const
    {
        Scalar mag = magnitude();
        if (mag < vSmallValue)
        {
            throw
                std::runtime_error
                (
                    "Error: Division by zero in Vector::normalized"
                );
        }

        Scalar inv = S(1.0) / mag;
        return Vector(x_ * inv, y_ * inv, z_ * inv);
    }

private:

    /// X, Y, Z components of the vector
    Scalar x_ = 0.0;
    Scalar y_ = 0.0;
    Scalar z_ = 0.0;
};

// Non-member methods

/**
 * @brief Computes dot product of two vectors
 * @param p1 First vector
 * @param p2 Second vector
 * @return Scalar dot product
 */
inline Scalar dot(const Vector& p1, const Vector& p2) noexcept
{
    return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}

/**
 * @brief Computes cross product of two vectors
 * @param p1 First vector
 * @param p2 Second vector
 * @return Cross product vector
 */
inline Vector cross(const Vector& p1, const Vector& p2) noexcept
{
    return
        Vector
        (
            p1.y() * p2.z() - p1.z() * p2.y(),
            p1.z() * p2.x() - p1.x() * p2.z(),
            p1.x() * p2.y() - p1.y() * p2.x()
        );
}

/**
 * @brief Scalar multiplication operator (scalar * vector)
 * @param scalar Scalar multiplier
 * @param p Vector to multiply
 * @return Scaled vector
 */
inline Vector operator*(Scalar scalar, const Vector& p) noexcept
{
    return p * scalar;
}

/**
 * @brief Stream output operator for Vector
 * @param os Output stream
 * @param p Vector to output
 * @return Reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Vector& p);
