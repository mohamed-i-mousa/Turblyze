#include "Vector.h"
#include <stdexcept>
#include <iostream>
#include <cmath>

Vector::Vector() : x(0.0), y(0.0), z(0.0) {}

Vector::Vector(Scalar x_val, Scalar y_val, Scalar z_val) 
    : x(x_val), y(y_val), z(z_val) {}

Vector Vector::operator+(const Vector& other) const 
{
    return Vector(x + other.x, y + other.y, z + other.z);
}

Vector Vector::operator-(const Vector& other) const 
{
    return Vector(x - other.x, y - other.y, z - other.z);
}

Vector Vector::operator*(Scalar scalar) const 
{
    return Vector(x * scalar, y * scalar, z * scalar);
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

    return Vector(x / scalar, y / scalar, z / scalar);
}

Vector& Vector::operator+=(const Vector& other) 
{
    x += other.x;
    y += other.y;
    z += other.z;

    return *this;
}

Vector& Vector::operator-=(const Vector& other) 
{
    x -= other.x;
    y -= other.y;
    z -= other.z;

    return *this;
}

Vector& Vector::operator*=(Scalar scalar) 
{
    x *= scalar;
    y *= scalar;
    z *= scalar;

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

    x /= scalar;
    y /= scalar;
    z /= scalar;

    return *this;
}

bool Vector::operator==(const Vector& other) const 
{
    return (std::abs(x - other.x) < equalityTolerance) 
        && (std::abs(y - other.y) < equalityTolerance)
        && (std::abs(z - other.z) < equalityTolerance);
}

bool Vector::operator!=(const Vector& other) const 
{
    return !(*this == other);
}

Scalar Vector::magnitudeSquared() const 
{
    return x * x + y * y + z * z;
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
        x /= mag;
        y /= mag;
        z /= mag;
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

Vector operator*(Scalar scalar, const Vector& p)
{
    return Vector(scalar * p.x, scalar * p.y, scalar * p.z);
}

std::ostream& operator<<(std::ostream& os, const Vector& p)
{
    os  << "(" << p.x << ", " << p.y << ", " << p.z << ")";
    return os;
}

Scalar dot(const Vector& p1, const Vector& p2)
{
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

Vector cross(const Vector& p1, const Vector& p2) 
{
    return  Vector
            (
                p1.y * p2.z - p1.z * p2.y,
                p1.z * p2.x - p1.x * p2.z,
                p1.x * p2.y - p1.y * p2.x
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
