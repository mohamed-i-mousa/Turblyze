/******************************************************************************
 * @file Vector.h
 * @brief 3D vector class for geometric and mathematical operations in CFD
 *
 * @details This header defines a 3D vector class that serves as the foundation
 * for all vector-based calculations in the CFD solver. The Vector class
 * provides essential mathematical operations required in the finite volume
 * discretization and mesh operations.
 *
 * @class Vector
 * - Components access and manipulation (x, y, z coordinates)
 * - Arithmetic operations (addition, subtraction, scalar multiplication)
 * - Vector operations (dot product, cross product, normalization)
 * - Comparison operators for sorting and equality testing
 * - Stream I/O operators for debugging
 *****************************************************************************/

#pragma once

#include <cmath>
#include <iosfwd>

#include "Scalar.h"
#include "ErrorHandler.h"


class Vector
{
public:

    /// Default constructor
    Vector() noexcept = default;

    /// Construct vector with specified components
    Vector(Scalar xValue, Scalar yValue, Scalar zValue) noexcept
    :
        x_(xValue),
        y_(yValue),
        z_(zValue)
    {}

// Setter methods

    /// Set X component
    void setX(Scalar xValue) noexcept { x_ = xValue; }

    /// Set Y component
    void setY(Scalar yValue) noexcept { y_ = yValue; }

    /// Set Z component
    void setZ(Scalar zValue) noexcept { z_ = zValue; }

// Accessor methods

    /// Get X component
    [[nodiscard]] Scalar x() const noexcept { return x_; }

    /// Get Y component
    [[nodiscard]] Scalar y() const noexcept { return y_; }

    /// Get Z component
    [[nodiscard]] Scalar z() const noexcept { return z_; }

// Operator methods

    /// Vector addition operator
    Vector operator+(const Vector& other) const noexcept
    {
        Vector result(*this);
        result += other;
        return result;
    }

    /// Vector subtraction operator
    Vector operator-(const Vector& other) const noexcept
    {
        Vector result(*this);
        result -= other;
        return result;
    }

    /// Scalar multiplication operator
    Vector operator*(Scalar scalar) const noexcept
    {
        Vector result(*this);
        result *= scalar;
        return result;
    }

    /// Scalar division operator
    Vector operator/(Scalar scalar) const noexcept
    {
        Vector result(*this);
        result /= scalar;
        return result;
    }

    /// Compound addition assignment operator
    Vector& operator+=(const Vector& other) noexcept
    {
        x_ += other.x_;
        y_ += other.y_;
        z_ += other.z_;

        return *this;
    }

    /// Compound subtraction assignment operator
    Vector& operator-=(const Vector& other) noexcept
    {
        x_ -= other.x_;
        y_ -= other.y_;
        z_ -= other.z_;

        return *this;
    }

    /// Compound multiplication assignment operator
    Vector& operator*=(Scalar scalar) noexcept
    {
        x_ *= scalar;
        y_ *= scalar;
        z_ *= scalar;

        return *this;
    }

    /// Compound division assignment operator
    Vector& operator/=(Scalar scalar) noexcept
    {
        if (std::abs(scalar) <= vSmallValue)
        {
            FatalError("Division by zero in Vector::operator/=");
        }

        const Scalar inverse = S(1.0) / scalar;
        x_ *= inverse;
        y_ *= inverse;
        z_ *= inverse;

        return *this;
    }

    /// Equality comparison operator
    bool operator==(const Vector& other) const noexcept
    {
        return (std::abs(x_ - other.x_) <= vSmallValue)
            && (std::abs(y_ - other.y_) <= vSmallValue)
            && (std::abs(z_ - other.z_) <= vSmallValue);
    }

private:

    /// x, y, z components of the vector
    Scalar x_ = S(0.0);
    Scalar y_ = S(0.0);
    Scalar z_ = S(0.0);
};

// Non-member methods

/// Compute dot product of two vectors
[[nodiscard]] inline Scalar dot(const Vector& p1, const Vector& p2) noexcept
{
    return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}

/// Compute cross product of two vectors
[[nodiscard]] inline Vector cross(const Vector& p1, const Vector& p2) noexcept
{
    return
        Vector
        (
            p1.y() * p2.z() - p1.z() * p2.y(),
            p1.z() * p2.x() - p1.x() * p2.z(),
            p1.x() * p2.y() - p1.y() * p2.x()
        );
}

/// Squared magnitude of a vector
[[nodiscard]] inline Scalar magnitudeSquared(const Vector& v) noexcept
{
    return v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
}

/// Magnitude of a vector
[[nodiscard]] inline Scalar magnitude(const Vector& v) noexcept
{
    return std::sqrt(magnitudeSquared(v));
}

/// Return a normalized copy of a vector
[[nodiscard]] inline Vector normalized(const Vector& v) noexcept
{
    const Scalar mag = magnitude(v);

    if (mag < vSmallValue)
    {
        FatalError("Division by zero in normalized(Vector)");
    }

    const Scalar inverse = S(1.0) / mag;
    return Vector(v.x() * inverse, v.y() * inverse, v.z() * inverse);
}

/// Scalar multiplication operator
inline Vector operator*(Scalar scalar, const Vector& p) noexcept
{
    return p * scalar;
}

/// Stream output operator for Vector
std::ostream& operator<<(std::ostream& os, const Vector& p);
