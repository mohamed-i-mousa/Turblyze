/******************************************************************************
 * @file BoundaryData.cpp
 * @brief Implementation of boundary condition data storage and retrieval
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "BoundaryData.h"

// Project headers
#include "ErrorHandler.h"

// ****************************** Setter Methods ******************************

void BoundaryData::setFixedValue(Scalar scalarValue) noexcept
{
    type_ = BCType::fixedValue;
    scalarValue_ = scalarValue;
}


void BoundaryData::setFixedGradient(Scalar scalarGradient) noexcept
{
    type_ = BCType::fixedGradient;
    scalarGradient_ = scalarGradient;
}


void BoundaryData::setZeroGradient() noexcept
{
    type_ = BCType::zeroGradient;
    scalarGradient_ = S(0.0);
}


void BoundaryData::setNoSlip() noexcept
{
    type_ = BCType::noSlip;
    scalarValue_ = S(0.0);
}


void BoundaryData::setWallFunctionType(BCType wallType) noexcept
{
    type_ = wallType;
}

// ***************************** Accessor Methods *****************************

Scalar BoundaryData::fixedScalarValue() const noexcept
{
    if (type_ == BCType::fixedValue || type_ == BCType::noSlip)
    {
        return scalarValue_;
    }

    FatalError
    (
        "Attempted to get fixed scalar value, but BC is not set to "
        "fixedValue or noSlip."
    );
}


Scalar BoundaryData::fixedScalarGradient() const noexcept
{
    if (type_ == BCType::fixedGradient)
    {
        return scalarGradient_;
    }

    FatalError
    (
        "Attempted to get fixed scalar gradient, "
        "but BC is not set to fixedGradient."
    );
}

// *************************** Non-Member Functions ***************************

NameRef bcTypeToString(BCType bctype) noexcept
{
    using enum BCType;
    switch (bctype)
    {
        case undefined:         return "undefined";
        case fixedValue:        return "fixedValue";
        case fixedGradient:     return "fixedGradient";
        case zeroGradient:      return "zeroGradient";
        case noSlip:            return "noSlip";
        case kWallFunction:     return "kWallFunction";
        case omegaWallFunction: return "omegaWallFunction";
        case nutWallFunction:   return "nutWallFunction";
    }

    FatalError("Corrupted BCType value");
}
