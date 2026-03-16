/******************************************************************************
 * @file BoundaryData.cpp
 * @brief Implementation of boundary condition data storage and retrieval
 *****************************************************************************/

#include "BoundaryData.hpp"
#include <stdexcept>


// ****************************** Setter Methods ******************************

void BoundaryData::setFixedValue(Scalar scalarValue) noexcept
{
    type_ = BCType::FIXED_VALUE;
    scalarValue_ = scalarValue;
    valueType_ = BCValueType::SCALAR;

    vectorValue_ = Vector();
    gradientType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setFixedValue(const Vector& vectorValue) noexcept
{
    type_ = BCType::FIXED_VALUE;
    vectorValue_ = vectorValue;
    valueType_ = BCValueType::VECTOR;

    scalarValue_ = S(0.0);
    gradientType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setFixedGradient(Scalar scalarGradient) noexcept
{
    type_ = BCType::FIXED_GRADIENT;
    scalarGradient_ = scalarGradient;
    gradientType_ = BCValueType::SCALAR;

    vectorValue_ = Vector();
    valueType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setFixedGradient(const Vector& vectorGradient) noexcept
{
    type_ = BCType::FIXED_GRADIENT;
    vectorGradient_ = vectorGradient;
    gradientType_ = BCValueType::VECTOR;

    scalarValue_ = S(0.0);
    valueType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setZeroGradient() noexcept
{
    type_ = BCType::ZERO_GRADIENT;
    scalarGradient_ = S(0.0);
    vectorGradient_ = Vector(S(0.0), S(0.0), S(0.0));

    valueType_ = BCValueType::UNDEFINED;
    gradientType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setNoSlip() noexcept
{
    type_ = BCType::NO_SLIP;
    vectorValue_ = Vector(S(0.0), S(0.0), S(0.0));
    valueType_ = BCValueType::VECTOR;

    scalarValue_ = S(0.0);
    gradientType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setKWallFunction() noexcept
{
    type_ = BCType::K_WALL_FUNCTION;
    scalarGradient_ = S(0.0);
    vectorGradient_ = Vector(S(0.0), S(0.0), S(0.0));

    valueType_ = BCValueType::UNDEFINED;
    gradientType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setOmegaWallFunction() noexcept
{
    type_ = BCType::OMEGA_WALL_FUNCTION;
    scalarGradient_ = S(0.0);
    vectorGradient_ = Vector(S(0.0), S(0.0), S(0.0));

    valueType_ = BCValueType::UNDEFINED;
    gradientType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setNutWallFunction() noexcept
{
    type_ = BCType::NUT_WALL_FUNCTION;
    scalarGradient_ = S(0.0);
    vectorGradient_ = Vector(S(0.0), S(0.0), S(0.0));

    valueType_ = BCValueType::UNDEFINED;
    gradientType_ = BCValueType::UNDEFINED;
}


// ***************************** Accessor Methods *****************************

Scalar BoundaryData::fixedScalarValue() const
{
    if (type_ == BCType::FIXED_VALUE && valueType_ == BCValueType::SCALAR)
    {
        return scalarValue_;
    }

    throw
        std::runtime_error
        (
            "Attempted to get fixed scalar value, "
            "but BC is not set to "
            "FIXED_VALUE with SCALAR type."
        );
}

const Vector& BoundaryData::fixedVectorValue() const
{
    if
    (
        (type_ == BCType::FIXED_VALUE || type_ == BCType::NO_SLIP)
     && valueType_ == BCValueType::VECTOR
    )
    {
        return vectorValue_;
    }

    throw
        std::runtime_error
        (
            "Attempted to get fixed vector value, "
            "but BC is not set to "
            "FIXED_VALUE/NO_SLIP with VECTOR type."
        );
}

Scalar BoundaryData::fixedScalarGradient() const
{
    if
    (
        type_ == BCType::FIXED_GRADIENT
     && gradientType_ == BCValueType::SCALAR
    )
    {
        return scalarGradient_;
    }

    throw
        std::runtime_error
        (
            "Attempted to get fixed scalar gradient, "
            "but BC is not set to "
            "FIXED_GRADIENT with SCALAR type."
        );
}

const Vector& BoundaryData::fixedVectorGradient() const
{
    if
    (
        type_ == BCType::FIXED_GRADIENT
     && gradientType_ == BCValueType::VECTOR
    )
    {
        return vectorGradient_;
    }

    throw
        std::runtime_error
        (
            "Attempted to get fixed vector gradient, "
            "but BC is not set to "
            "FIXED_GRADIENT with VECTOR type."
        );
}