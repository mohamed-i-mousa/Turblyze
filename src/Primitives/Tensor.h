/******************************************************************************
 * @file Tensor.h
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

// ********************************** Headers *********************************

// Standard library headers
#include <cmath>
#include <iosfwd>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "ErrorHandler.h"
#include "Integer.h"

// ******************************* class Tensor *******************************

class Tensor
{
public:

    /// Default constructor (zero-initialised components)
    Tensor() noexcept = default;

    /// Constructs tensor with specified components (row-major)
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

// ****************************** Setter Methods ******************************

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

// ***************************** Accessor Methods *****************************

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

    /// Get row i of the tensor
    [[nodiscard]] Vector row(Index i) const noexcept
    {
        switch (i)
        {
            case 0: return Vector{xx_, xy_, xz_};
            case 1: return Vector{yx_, yy_, yz_};
            case 2: return Vector{zx_, zy_, zz_};
        }

        FatalError("Tensor::row index out of range");
        return Vector{};
    }

    /// Get column j of the tensor
    [[nodiscard]] Vector col(Index j) const noexcept
    {
        switch (j)
        {
            case 0: return Vector{xx_, yx_, zx_};
            case 1: return Vector{xy_, yy_, zy_};
            case 2: return Vector{xz_, yz_, zz_};
        }

        FatalError("Tensor::col index out of range");
        return Vector{};
    }

// ***************************** Operator Methods *****************************

    /// Tensor addition operator
    Tensor operator+(const Tensor& other) const noexcept
    {
        Tensor result(*this);
        result += other;
        return result;
    }

    /// Tensor subtraction operator
    Tensor operator-(const Tensor& other) const noexcept
    {
        Tensor result(*this);
        result -= other;
        return result;
    }

    /// Scalar multiplication operator
    Tensor operator*(Scalar scalar) const noexcept
    {
        Tensor result(*this);
        result *= scalar;
        return result;
    }

    /// Scalar division operator
    Tensor operator/(Scalar scalar) const noexcept
    {
        Tensor result(*this);
        result /= scalar;
        return result;
    }

    /// Compound addition assignment operator
    Tensor& operator+=(const Tensor& other) noexcept
    {
        xx_ += other.xx_; xy_ += other.xy_; xz_ += other.xz_;
        yx_ += other.yx_; yy_ += other.yy_; yz_ += other.yz_;
        zx_ += other.zx_; zy_ += other.zy_; zz_ += other.zz_;

        return *this;
    }

    /// Compound subtraction assignment operator
    Tensor& operator-=(const Tensor& other) noexcept
    {
        xx_ -= other.xx_; xy_ -= other.xy_; xz_ -= other.xz_;
        yx_ -= other.yx_; yy_ -= other.yy_; yz_ -= other.yz_;
        zx_ -= other.zx_; zy_ -= other.zy_; zz_ -= other.zz_;

        return *this;
    }

    /// Compound multiplication assignment operator
    Tensor& operator*=(Scalar scalar) noexcept
    {
        xx_ *= scalar; xy_ *= scalar; xz_ *= scalar;
        yx_ *= scalar; yy_ *= scalar; yz_ *= scalar;
        zx_ *= scalar; zy_ *= scalar; zz_ *= scalar;

        return *this;
    }

    /// Compound division assignment operator
    Tensor& operator/=(Scalar scalar) noexcept
    {
        if (std::abs(scalar) <= vSmallValue)
        {
            FatalError("Division by zero in Tensor::operator/=");
        }

        const Scalar inverse = S(1.0) / scalar;
        xx_ *= inverse; xy_ *= inverse; xz_ *= inverse;
        yx_ *= inverse; yy_ *= inverse; yz_ *= inverse;
        zx_ *= inverse; zy_ *= inverse; zz_ *= inverse;

        return *this;
    }

// ****************************** Tensor Algebra ******************************

    /// Transpose of this tensor
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

    /// Symmetric part of this tensor: 0.5 * (T + T^T)
    [[nodiscard]] Tensor symm() const noexcept
    {
        const Scalar sxy = S(0.5) * (xy_ + yx_);
        const Scalar sxz = S(0.5) * (xz_ + zx_);
        const Scalar syz = S(0.5) * (yz_ + zy_);

        return
            Tensor
            (
                xx_, sxy, sxz,
                sxy, yy_, syz,
                sxz, syz, zz_
            );
    }

    /// Antisymmetric (skew) part of this tensor: 0.5 * (T - T^T)
    [[nodiscard]] Tensor skew() const noexcept
    {
        const Scalar axy = S(0.5) * (xy_ - yx_);
        const Scalar axz = S(0.5) * (xz_ - zx_);
        const Scalar ayz = S(0.5) * (yz_ - zy_);

        return
            Tensor
            (
                S(0.0),  axy,     axz,
                -axy,    S(0.0),  ayz,
                -axz,    -ayz,    S(0.0)
            );
    }

    /// Trace of this tensor
    [[nodiscard]] Scalar trace() const noexcept
    {
        return xx_ + yy_ + zz_;
    }

    /// Frobenius squared magnitude: T:T = T_ij T_ij
    [[nodiscard]] Scalar magnitudeSquared() const noexcept
    {
        return 
            xx_ * xx_ + xy_ * xy_ + xz_ * xz_
          + yx_ * yx_ + yy_ * yy_ + yz_ * yz_
          + zx_ * zx_ + zy_ * zy_ + zz_ * zz_;
    }

// ****************************** Private Members *****************************

private:

    /// Row-major components (first index = row, second = column)
    Scalar xx_ = S(0.0); Scalar xy_ = S(0.0); Scalar xz_ = S(0.0);
    Scalar yx_ = S(0.0); Scalar yy_ = S(0.0); Scalar yz_ = S(0.0);
    Scalar zx_ = S(0.0); Scalar zy_ = S(0.0); Scalar zz_ = S(0.0);
};

// *************************** Non-Member Functions ***************************

/// Construct a tensor from three row vectors
[[nodiscard]] inline Tensor tensorFromRows
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

/// Double-dot product of two tensors: A:B = A_ij B_ij
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

/// Outer product of two vectors: T_ij = a_i b_j
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

/// Scalar multiplication operator (scalar * tensor)
inline Tensor operator*(Scalar scalar, const Tensor& T) noexcept
{
    return T * scalar;
}

/// Stream output operator for Tensor
std::ostream& operator<<(std::ostream& os, const Tensor& T);
