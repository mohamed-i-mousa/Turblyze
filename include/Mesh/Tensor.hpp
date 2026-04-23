/******************************************************************************
 * @file Tensor.hpp
 * @brief 3x3 tensor class for geometric and mathematical operations
 *
 * @details This header defines a 3x3 tensor class used for velocity
 * gradients, strain-rate and rotation tensors, Reynolds-stress tensors,
 * and any other second-order tensor fields in the solver.
 *
 * @class Tensor
 * - Component access and manipulation (row-major xx..zz)
 * - Arithmetic operations (addition, subtraction, scalar multiplication)
 * - Tensor operations (transpose, symmetric/antisymmetric part, trace)
 * - Double-dot product and Frobenius magnitude
 * - Comparison operators and stream I/O for debugging
 *****************************************************************************/

#pragma once

#include <cmath>
#include <cstddef>
#include <iosfwd>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "ErrorHandler.hpp"


class Tensor
{
public:

    /// Default constructor (zero-initialised components)
    Tensor() noexcept = default;

    /**
     * @brief Constructs tensor with specified components (row-major)
     * @param xx Component (0,0)
     * @param xy Component (0,1)
     * @param xz Component (0,2)
     * @param yx Component (1,0)
     * @param yy Component (1,1)
     * @param yz Component (1,2)
     * @param zx Component (2,0)
     * @param zy Component (2,1)
     * @param zz Component (2,2)
     */
    Tensor
    (
        Scalar xx, Scalar xy, Scalar xz,
        Scalar yx, Scalar yy, Scalar yz,
        Scalar zx, Scalar zy, Scalar zz
    ) noexcept
    :
        xx_(xx), xy_(xy), xz_(xz),
        yx_(yx), yy_(yy), yz_(yz),
        zx_(zx), zy_(zy), zz_(zz)
    {}

// Factory methods

    /**
     * @brief Constructs tensor from three row vectors
     * @param row0 First row (components at i=0)
     * @param row1 Second row (components at i=1)
     * @param row2 Third row (components at i=2)
     * @return Tensor with the given rows
     */
    [[nodiscard]] static Tensor fromRows
    (
        const Vector& row0,
        const Vector& row1,
        const Vector& row2
    ) noexcept
    {
        return Tensor
        (
            row0.x(), row0.y(), row0.z(),
            row1.x(), row1.y(), row1.z(),
            row2.x(), row2.y(), row2.z()
        );
    }

// Setter methods

    /// Set component (0,0)
    void setXX(Scalar value) noexcept { xx_ = value; }

    /// Set component (0,1)
    void setXY(Scalar value) noexcept { xy_ = value; }

    /// Set component (0,2)
    void setXZ(Scalar value) noexcept { xz_ = value; }

    /// Set component (1,0)
    void setYX(Scalar value) noexcept { yx_ = value; }

    /// Set component (1,1)
    void setYY(Scalar value) noexcept { yy_ = value; }

    /// Set component (1,2)
    void setYZ(Scalar value) noexcept { yz_ = value; }

    /// Set component (2,0)
    void setZX(Scalar value) noexcept { zx_ = value; }

    /// Set component (2,1)
    void setZY(Scalar value) noexcept { zy_ = value; }

    /// Set component (2,2)
    void setZZ(Scalar value) noexcept { zz_ = value; }

// Accessor methods

    /// Get component (0,0)
    [[nodiscard]] Scalar xx() const noexcept { return xx_; }

    /// Get component (0,1)
    [[nodiscard]] Scalar xy() const noexcept { return xy_; }

    /// Get component (0,2)
    [[nodiscard]] Scalar xz() const noexcept { return xz_; }

    /// Get component (1,0)
    [[nodiscard]] Scalar yx() const noexcept { return yx_; }

    /// Get component (1,1)
    [[nodiscard]] Scalar yy() const noexcept { return yy_; }

    /// Get component (1,2)
    [[nodiscard]] Scalar yz() const noexcept { return yz_; }

    /// Get component (2,0)
    [[nodiscard]] Scalar zx() const noexcept { return zx_; }

    /// Get component (2,1)
    [[nodiscard]] Scalar zy() const noexcept { return zy_; }

    /// Get component (2,2)
    [[nodiscard]] Scalar zz() const noexcept { return zz_; }

    /**
     * @brief Get row i of the tensor
     * @param i Row index (0, 1, or 2)
     * @return Vector containing the three components of row i
     * @note Terminates the program if i is out of range
     */
    [[nodiscard]] Vector row(size_t i) const
    {
        switch (i)
        {
            case 0: return Vector(xx_, xy_, xz_);
            case 1: return Vector(yx_, yy_, yz_);
            case 2: return Vector(zx_, zy_, zz_);
        }

        FatalError("Tensor::row index out of range");
        return Vector{};
    }

    /**
     * @brief Get column j of the tensor
     * @param j Column index (0, 1, or 2)
     * @return Vector containing the three components of column j
     * @note Terminates the program if j is out of range
     */
    [[nodiscard]] Vector col(size_t j) const
    {
        switch (j)
        {
            case 0: return Vector(xx_, yx_, zx_);
            case 1: return Vector(xy_, yy_, zy_);
            case 2: return Vector(xz_, yz_, zz_);
        }

        FatalError("Tensor::col index out of range");
        return Vector{};
    }

// Operator methods

    /**
     * @brief Tensor addition operator
     * @param other Tensor to add
     * @return Sum of the two tensors
     */
    Tensor operator+(const Tensor& other) const noexcept
    {
        Tensor result(*this);
        result += other;
        return result;
    }

    /**
     * @brief Tensor subtraction operator
     * @param other Tensor to subtract
     * @return Difference of the two tensors
     */
    Tensor operator-(const Tensor& other) const noexcept
    {
        Tensor result(*this);
        result -= other;
        return result;
    }

    /**
     * @brief Scalar multiplication operator
     * @param scalar Scalar to multiply by
     * @return Tensor scaled by scalar
     */
    Tensor operator*(Scalar scalar) const noexcept
    {
        Tensor result(*this);
        result *= scalar;
        return result;
    }

    /**
     * @brief Scalar division operator
     * @param scalar Scalar to divide by
     * @return Tensor divided by scalar
     * @note Terminates the program if scalar is near zero
     */
    Tensor operator/(Scalar scalar) const
    {
        Tensor result(*this);
        result /= scalar;
        return result;
    }

    /**
     * @brief Compound addition assignment operator
     * @param other Tensor to add
     * @return Reference to this tensor
     */
    Tensor& operator+=(const Tensor& other) noexcept
    {
        xx_ += other.xx_; xy_ += other.xy_; xz_ += other.xz_;
        yx_ += other.yx_; yy_ += other.yy_; yz_ += other.yz_;
        zx_ += other.zx_; zy_ += other.zy_; zz_ += other.zz_;

        return *this;
    }

    /**
     * @brief Compound subtraction assignment operator
     * @param other Tensor to subtract
     * @return Reference to this tensor
     */
    Tensor& operator-=(const Tensor& other) noexcept
    {
        xx_ -= other.xx_; xy_ -= other.xy_; xz_ -= other.xz_;
        yx_ -= other.yx_; yy_ -= other.yy_; yz_ -= other.yz_;
        zx_ -= other.zx_; zy_ -= other.zy_; zz_ -= other.zz_;

        return *this;
    }

    /**
     * @brief Compound multiplication assignment operator
     * @param scalar Scalar to multiply by
     * @return Reference to this tensor
     */
    Tensor& operator*=(Scalar scalar) noexcept
    {
        xx_ *= scalar; xy_ *= scalar; xz_ *= scalar;
        yx_ *= scalar; yy_ *= scalar; yz_ *= scalar;
        zx_ *= scalar; zy_ *= scalar; zz_ *= scalar;

        return *this;
    }

    /**
     * @brief Compound division assignment operator
     * @param scalar Scalar to divide by
     * @return Reference to this tensor
     * @note Terminates the program if scalar is near zero
     * @note Using inverse is a micro-optimization to reduce the number
     *       of divisions, which are more expensive than multiplications
     */
    Tensor& operator/=(Scalar scalar)
    {
        if (std::abs(scalar) <= vSmallValue)
        {
            FatalError("Division by zero in Tensor::operator/=");
        }

        Scalar inverse = S(1.0) / scalar;
        xx_ *= inverse; xy_ *= inverse; xz_ *= inverse;
        yx_ *= inverse; yy_ *= inverse; yz_ *= inverse;
        zx_ *= inverse; zy_ *= inverse; zz_ *= inverse;

        return *this;
    }

// Tensor algebra

    /**
     * @brief Transpose of this tensor
     * @return Tensor with components (i,j) and (j,i) swapped
     */
    [[nodiscard]] Tensor transpose() const noexcept
    {
        return
            Tensor
            (
                xx_, yx_, zx_,
                xy_, yy_, zy_,
                xz_, yz_, zz_
            );
    }

    /**
     * @brief Symmetric part of this tensor: 0.5 * (T + T^T)
     * @return Symmetric tensor
     */
    [[nodiscard]] Tensor symm() const noexcept
    {
        Scalar sxy = S(0.5) * (xy_ + yx_);
        Scalar sxz = S(0.5) * (xz_ + zx_);
        Scalar syz = S(0.5) * (yz_ + zy_);

        return
            Tensor
            (
                xx_, sxy, sxz,
                sxy, yy_, syz,
                sxz, syz, zz_
            );
    }

    /**
     * @brief Antisymmetric (skew) part of this tensor: 0.5 * (T - T^T)
     * @return Skew-symmetric tensor (diagonals are zero)
     */
    [[nodiscard]] Tensor skew() const noexcept
    {
        Scalar axy = S(0.5) * (xy_ - yx_);
        Scalar axz = S(0.5) * (xz_ - zx_);
        Scalar ayz = S(0.5) * (yz_ - zy_);

        return
            Tensor
            (
                S(0.0),  axy,     axz,
                -axy,    S(0.0),  ayz,
                -axz,    -ayz,    S(0.0)
            );
    }

    /**
     * @brief Trace of this tensor
     * @return Sum of diagonal components (xx + yy + zz)
     */
    [[nodiscard]] Scalar trace() const noexcept
    {
        return xx_ + yy_ + zz_;
    }

    /**
     * @brief Frobenius squared magnitude: T:T = T_ij T_ij
     * @return Sum of squared components
     */
    [[nodiscard]] Scalar magnitudeSquared() const noexcept
    {
        return 
            xx_ * xx_ + xy_ * xy_ + xz_ * xz_
          + yx_ * yx_ + yy_ * yy_ + yz_ * yz_
          + zx_ * zx_ + zy_ * zy_ + zz_ * zz_;
    }

private:

    /// Row-major components (first index = row, second = column)
    Scalar xx_ = 0.0; Scalar xy_ = 0.0; Scalar xz_ = 0.0;
    Scalar yx_ = 0.0; Scalar yy_ = 0.0; Scalar yz_ = 0.0;
    Scalar zx_ = 0.0; Scalar zy_ = 0.0; Scalar zz_ = 0.0;
};

// Non-member methods

/**
 * @brief Double-dot product of two tensors: A:B = A_ij B_ij
 * @param A First tensor
 * @param B Second tensor
 * @return Scalar double-dot product
 */
[[nodiscard]] inline Scalar doubleDot
(
    const Tensor& A,
    const Tensor& B
) noexcept
{
    return 
        A.xx() * B.xx() + A.xy() * B.xy() + A.xz() * B.xz()
      + A.yx() * B.yx() + A.yy() * B.yy() + A.yz() * B.yz()
      + A.zx() * B.zx() + A.zy() * B.zy() + A.zz() * B.zz();
}

/**
 * @brief Outer product of two vectors: T_ij = a_i b_j
 * @param a First vector
 * @param b Second vector
 * @return Tensor outer product
 */
[[nodiscard]] inline Tensor outer
(
    const Vector& a,
    const Vector& b
) noexcept
{
    return 
        Tensor
        (
            a.x() * b.x(), a.x() * b.y(), a.x() * b.z(),
            a.y() * b.x(), a.y() * b.y(), a.y() * b.z(),
            a.z() * b.x(), a.z() * b.y(), a.z() * b.z()
        );
}

/**
 * @brief Scalar multiplication operator (scalar * tensor)
 * @param scalar Scalar multiplier
 * @param T Tensor to multiply
 * @return Scaled tensor
 */
inline Tensor operator*(Scalar scalar, const Tensor& T) noexcept
{
    return T * scalar;
}

/**
 * @brief Stream output operator for Tensor
 * @param os Output stream
 * @param T Tensor to output
 * @return Reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Tensor& T);
