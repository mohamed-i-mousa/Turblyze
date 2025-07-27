#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <limits>

#include "Scalar.h"

struct Vector {
    // ----- Members ----- //
    Scalar x, y, z;

    // ----- Constructors ----- //

    // Default constructor
    Vector() : x(0.0), y(0.0), z(0.0) {}

    // Parametized constructor
    Vector(Scalar x_val, Scalar y_val, Scalar z_val) : x(x_val), y(y_val), z(z_val) {}

    // ----- Operators ----- //

    // Addition of two vectors
    Vector operator+(const Vector& other) const {
        return Vector(x + other.x, y + other.y, z + other.z);
    }

    // Subtraction of two vectors
    Vector operator-(const Vector& other) const {
        return Vector(x - other.x, y - other.y, z - other.z);
    }

    // Multiplication of a vector with a scalar
    Vector operator*(Scalar scalar) const {
        return Vector(x * scalar, y * scalar, z * scalar);
    }

    // Division of a vector by a scalar
    Vector operator/(Scalar scalar) const {
        if (std::abs(scalar) < DIVISION_TOLERANCE) {
            throw std::runtime_error("Error: Division by zero in Vector::operator/");
        }
        return Vector(x / scalar, y / scalar, z / scalar);
    }

    // Compound assignment operators
    Vector& operator+=(const Vector& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Vector& operator-=(const Vector& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    Vector& operator*=(Scalar scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    Vector& operator/=(Scalar scalar) {
        if (std::abs(scalar) < DIVISION_TOLERANCE) {
            throw std::runtime_error("Error: Division by zero in Vector::operator/=");
        }
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }

    // Equality comparison: Vector1 == Vector2
    // Uses a small tolerance for floating-point comparisons
    bool operator==(const Vector& other) const {
        return (std::abs(x - other.x) < EQUALITY_TOLERANCE) &&
               (std::abs(y - other.y) < EQUALITY_TOLERANCE) &&
               (std::abs(z - other.z) < EQUALITY_TOLERANCE);
    }

    // Inequality comparison: Vector1 != Vector2
    bool operator!=(const Vector& other) const {
        return !(*this == other);
    }

    // ----- Methods ----- //

    // Calculate the squared magnitude
    Scalar magnitudeSquared() const {
        return x * x + y * y + z * z;
    }

    // Calculate the magnitude
    Scalar magnitude() const {
        return std::sqrt(magnitudeSquared());
    }

    // Normalize the vector
    // Modifies the vector and returns a reference to itself
    Vector& normalize() {
        Scalar mag = magnitude();
        if (std::abs(mag) > DIVISION_TOLERANCE) {
            x /= mag;
            y /= mag;
            z /= mag;
        }
        else {
            throw std::runtime_error("Error: Division by zero in Vector::normalize");
        }

        return *this;
    }

    // Get a normalized copy
    Vector normalized() const {
        Vector result = *this;
        return result.normalize();
    }
};

    // ----- Operator Overloads (Non-Member Methods) ----- //

    // Scalar Multiplication: scalar * P (allows Scalar * Vector syntax)
    inline Vector operator*(Scalar scalar, const Vector& p) {
        return Vector(scalar * p.x, scalar * p.y, scalar * p.z);
    }

    // Overload for std::ostream to print Vector objects: cout << P
    inline std::ostream& operator<<(std::ostream& os, const Vector& p) {
        os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
        return os;
    }

    // ----- Common Geometric Methods (Non-Member) ----- //

    // Dot product of two vectors (vectors from origin)
    inline Scalar dot(const Vector& p1, const Vector& p2) {
        return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
    }

    // Cross product of two vectors (vectors from origin)
    inline Vector cross(const Vector& p1, const Vector& p2) {
        return Vector(
        p1.y * p2.z - p1.z * p2.y,
        p1.z * p2.x - p1.x * p2.z,
        p1.x * p2.y - p1.y * p2.x
        );
    }

    // Distance between two vectors
    inline Scalar distance(const Vector& p1, const Vector& p2) {
        return (p2 - p1).magnitude();
    }

    // Squared distance between two vectors
    inline Scalar distanceSquared(const Vector& p1, const Vector& p2) {
        return (p2 - p1).magnitudeSquared();
    }

#endif
