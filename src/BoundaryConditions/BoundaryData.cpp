/******************************************************************************
 * @file BoundaryData.cpp
 * @brief Implementation of boundary condition data storage and retrieval
 *****************************************************************************/

#include "BoundaryData.h"

#include "ErrorHandler.h"


// ****************************** Setter Methods ******************************

void BoundaryData::setFixedValue(Scalar scalarValue) noexcept
{
    type_ = BCType::FIXED_VALUE;
    scalarValue_ = scalarValue;
    valueType_ = BCValueType::SCALAR;

    vectorValue_ = Vector{};
}

void BoundaryData::setFixedValue(const Vector& vectorValue) noexcept
{
    type_ = BCType::FIXED_VALUE;
    vectorValue_ = vectorValue;
    valueType_ = BCValueType::VECTOR;

    scalarValue_ = S(0.0);
}

void BoundaryData::setFixedGradient(Scalar scalarGradient) noexcept
{
    type_ = BCType::FIXED_GRADIENT;
    scalarGradient_ = scalarGradient;

    vectorValue_ = Vector{};
    valueType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setZeroGradient() noexcept
{
    type_ = BCType::ZERO_GRADIENT;
    scalarGradient_ = S(0.0);

    valueType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setNoSlip() noexcept
{
    type_ = BCType::NO_SLIP;
    vectorValue_ = Vector(S(0.0), S(0.0), S(0.0));
    valueType_ = BCValueType::VECTOR;

    scalarValue_ = S(0.0);
}

void BoundaryData::setWallFunctionType(BCType wallType) noexcept
{
    type_ = wallType;
    valueType_ = BCValueType::UNDEFINED;
}


// ***************************** Accessor Methods *****************************

Scalar BoundaryData::fixedScalarValue() const
{
    if (type_ == BCType::FIXED_VALUE && valueType_ == BCValueType::SCALAR)
    {
        return scalarValue_;
    }

    FatalError
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

    FatalError
    (
        "Attempted to get fixed vector value, "
        "but BC is not set to "
        "FIXED_VALUE/NO_SLIP with VECTOR type."
    );
}

Scalar BoundaryData::fixedScalarGradient() const
{
    if (type_ == BCType::FIXED_GRADIENT)
    {
        return scalarGradient_;
    }

    FatalError
    (
        "Attempted to get fixed scalar gradient, "
        "but BC is not set to "
        "FIXED_GRADIENT with SCALAR type."
    );
}
