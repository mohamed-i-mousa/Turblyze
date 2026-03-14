/******************************************************************************
 * @file Vector.cpp
 * @brief Implementation of 3D vector operations and utilities
 *****************************************************************************/

#include "Vector.hpp"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>

Vector::Vector(Scalar xValue, Scalar yValue, Scalar zValue) noexcept
    : x_(xValue), y_(yValue), z_(zValue) {}

// ************************** Operator Overloads **************************

Vector Vector::operator+(const Vector& other) const noexcept
{
    return Vector(x_ + other.x_, y_ + other.y_, z_ + other.z_);
}

Vector Vector::operator-(const Vector& other) const noexcept
{
    return Vector(x_ - other.x_, y_ - other.y_, z_ - other.z_);
}

Vector Vector::operator*(Scalar scalar) const noexcept
{
    return Vector(x_ * scalar, y_ * scalar, z_ * scalar);
}

Vector Vector::operator/(Scalar scalar) const
{
    if (std::abs(scalar) < vSmallValue)
    {
        throw
            std::runtime_error
            (
                "Error: Division by zero in Vector::operator/"
            );
    }

    return Vector(x_ / scalar, y_ / scalar, z_ / scalar);
}

Vector& Vector::operator+=(const Vector& other) noexcept
{
    x_ += other.x_;
    y_ += other.y_;
    z_ += other.z_;

    return *this;
}

Vector& Vector::operator-=(const Vector& other) noexcept
{
    x_ -= other.x_;
    y_ -= other.y_;
    z_ -= other.z_;

    return *this;
}

Vector& Vector::operator*=(Scalar scalar) noexcept
{
    x_ *= scalar;
    y_ *= scalar;
    z_ *= scalar;

    return *this;
}

Vector& Vector::operator/=(Scalar scalar)
{
    if (std::abs(scalar) < vSmallValue)
    {
        throw
            std::runtime_error
            (
                "Error: Division by zero in Vector::operator/="
            );
    }

    x_ /= scalar;
    y_ /= scalar;
    z_ /= scalar;

    return *this;
}

bool Vector::operator==(const Vector& other) const noexcept
{
    return (std::abs(x_ - other.x_) < smallValue)
        && (std::abs(y_ - other.y_) < smallValue)
        && (std::abs(z_ - other.z_) < smallValue);
}

bool Vector::operator!=(const Vector& other) const noexcept
{
    return !(*this == other);
}

Scalar Vector::magnitudeSquared() const noexcept
{
    return x_ * x_ + y_ * y_ + z_ * z_;
}

Scalar Vector::magnitude() const noexcept
{
    return std::sqrt(magnitudeSquared());
}

Vector Vector::normalized() const
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

    return Vector(x_ / mag, y_ / mag, z_ / mag);
}

// ************************* Non-Member Functions *************************

Vector operator*(Scalar scalar, const Vector& p) noexcept
{
    return p * scalar;
}

std::ostream& operator<<(std::ostream& os, const Vector& p)
{
    std::ios_base::fmtflags flags = os.flags();
    int prec = os.precision();

    os  << std::fixed << std::setprecision(6);

    os  << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";

    os.flags(flags);
    os.precision(prec);

    return os;
}

Scalar dot(const Vector& p1, const Vector& p2) noexcept
{
    return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}

Vector cross(const Vector& p1, const Vector& p2) noexcept
{
    return
        Vector
        (
            p1.y() * p2.z() - p1.z() * p2.y(),
            p1.z() * p2.x() - p1.x() * p2.z(),
            p1.x() * p2.y() - p1.y() * p2.x()
        );
}
