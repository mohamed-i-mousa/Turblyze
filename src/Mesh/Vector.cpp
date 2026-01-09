/******************************************************************************
 * @file Vector.cpp
 * @brief Implementation of 3D vector operations and utilities
 *****************************************************************************/

#include "Vector.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>

Vector::Vector() : x_(0.0), y_(0.0), z_(0.0) {}

Vector::Vector(Scalar x_value, Scalar y_value, Scalar z_value) 
    : x_(x_value), y_(y_value), z_(z_value) {}

// Public methods

Vector Vector::operator+(const Vector& other) const 
{
    return Vector(x_ + other.x_, y_ + other.y_, z_ + other.z_);
}

Vector Vector::operator-(const Vector& other) const 
{
    return Vector(x_ - other.x_, y_ - other.y_, z_ - other.z_);
}

Vector Vector::operator*(Scalar scalar) const 
{
    return Vector(x_ * scalar, y_ * scalar, z_ * scalar);
}

Vector Vector::operator/(Scalar scalar) const 
{
    if (std::abs(scalar) < smallValue) 
    {
        throw std::runtime_error
            (
                "Error: Division by zero in Vector::operator/"
            );
    }

    return Vector(x_ / scalar, y_ / scalar, z_ / scalar);
}

Vector& Vector::operator+=(const Vector& other) 
{
    x_ += other.x_;
    y_ += other.y_;
    z_ += other.z_;

    return *this;
}

Vector& Vector::operator-=(const Vector& other) 
{
    x_ -= other.x_;
    y_ -= other.y_;
    z_ -= other.z_;

    return *this;
}

Vector& Vector::operator*=(Scalar scalar) 
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
        throw std::runtime_error
            (
                "Error: Division by zero in Vector::operator/="
            );
    }

    x_ /= scalar;
    y_ /= scalar;
    z_ /= scalar;

    return *this;
}

bool Vector::operator==(const Vector& other) const 
{
    return (std::abs(x_ - other.x_) < smallValue) 
        && (std::abs(y_ - other.y_) < smallValue)
        && (std::abs(z_ - other.z_) < smallValue);
}

bool Vector::operator!=(const Vector& other) const 
{
    return !(*this == other);
}

Scalar Vector::magnitudeSquared() const 
{
    return x_ * x_ + y_ * y_ + z_ * z_;
}

Scalar Vector::magnitude() const 
{
    return std::sqrt(magnitudeSquared());
}

Vector& Vector::normalize()
{
    Scalar mag = magnitude();
    if (std::abs(mag) > smallValue)
    {
        x_ /= mag;
        y_ /= mag;
        z_ /= mag;
    }
    else
    {
        throw std::runtime_error(
            "Error: Division by zero in Vector::normalize");
    }

    return *this;
}

Vector Vector::normalized() const
{
    Vector result = *this;

    return result.normalize();
}

// Non-member methods 

Vector operator*(Scalar scalar, const Vector& p)
{
    return Vector(scalar * p.x(), scalar * p.y(), scalar * p.z());
}

std::ostream& operator<<(std::ostream& os, const Vector& p)
{
    os  << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";

    return os;
}

Scalar dot(const Vector& p1, const Vector& p2)
{
    return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}

Vector cross(const Vector& p1, const Vector& p2) 
{
    return  Vector
            (
                p1.y() * p2.z() - p1.z() * p2.y(),
                p1.z() * p2.x() - p1.x() * p2.z(),
                p1.x() * p2.y() - p1.y() * p2.x()
            );
}

Scalar distance(const Vector& p1, const Vector& p2)
{
    return (p2 - p1).magnitude();
}

Scalar distanceSquared(const Vector& p1, const Vector& p2)
{
    return (p2 - p1).magnitudeSquared();
}
