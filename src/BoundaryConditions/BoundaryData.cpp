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
}


void BoundaryData::setFixedGradient(Scalar scalarGradient) noexcept
{
    type_ = BCType::FIXED_GRADIENT;
    scalarGradient_ = scalarGradient;
}


void BoundaryData::setZeroGradient() noexcept
{
    type_ = BCType::ZERO_GRADIENT;
    scalarGradient_ = S(0.0);
}


void BoundaryData::setNoSlip() noexcept
{
    type_ = BCType::NO_SLIP;
    scalarValue_ = S(0.0);
}


void BoundaryData::setWallFunctionType(BCType wallType) noexcept
{
    type_ = wallType;
}

// ***************************** Accessor Methods *****************************

Scalar BoundaryData::fixedScalarValue() const noexcept
{
    if (type_ == BCType::FIXED_VALUE)
    {
        return scalarValue_;
    }

    FatalError
    (
        "Attempted to get fixed scalar value, "
        "but BC is not set to FIXED_VALUE."
    );
}


Scalar BoundaryData::fixedScalarGradient() const noexcept
{
    if (type_ == BCType::FIXED_GRADIENT)
    {
        return scalarGradient_;
    }

    FatalError
    (
        "Attempted to get fixed scalar gradient, "
        "but BC is not set to FIXED_GRADIENT."
    );
}
