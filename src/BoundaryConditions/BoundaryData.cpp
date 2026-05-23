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

// *************************** Non-Member Functions ***************************

std::string_view bcTypeToString(BCType bctype) noexcept
{
    using enum BCType;
    switch (bctype)
    {
        case UNDEFINED:           return "UNDEFINED";
        case FIXED_VALUE:         return "FIXED_VALUE";
        case FIXED_GRADIENT:      return "FIXED_GRADIENT";
        case ZERO_GRADIENT:       return "ZERO_GRADIENT";
        case NO_SLIP:             return "NO_SLIP";
        case K_WALL_FUNCTION:     return "K_WALL_FUNCTION";
        case OMEGA_WALL_FUNCTION: return "OMEGA_WALL_FUNCTION";
        case NUT_WALL_FUNCTION:   return "NUT_WALL_FUNCTION";
    }

    FatalError("Corrupted BCType value");
}
