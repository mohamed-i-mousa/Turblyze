/******************************************************************************
 * @file Constraint.cpp
 * @brief Implementation of constraint handling for CFD fields
 *****************************************************************************/

#include "Constraint.hpp"

#include <cmath>

#include "ErrorHandler.hpp"

Constraint::Constraint
(
    VectorField& velocityField,
    ScalarField& pressureField
) noexcept
:
    U_(velocityField),
    p_(pressureField)
{}


// ****************************** Setter Methods ******************************

void Constraint::setPressureConstraints
(
    Scalar minPressure,
    Scalar maxPressure
)
{
    if (minPressure < maxPressure)
    {
        minPressure_ = minPressure;
        maxPressure_ = maxPressure;
    }
    else if (minPressure > maxPressure)
    {
        minPressure_ = maxPressure;
        maxPressure_ = minPressure;
    }
    else
    {
        FatalError("Wrong Pressure Constraints input");
    }
}

void Constraint::enableConstraints
(
    bool enableVel,
    bool enablePress
) noexcept
{
    enableVelocityConstraints_ = enableVel;
    enablePressureConstraints_ = enablePress;
}

size_t Constraint::applyVelocityConstraints() noexcept
{
    if (!enableVelocityConstraints_) return 0;

    size_t constraintApplied = 0;
    Scalar maxMagSq = maxVelocityMagnitude_ * maxVelocityMagnitude_;

    for (size_t cellIdx = 0; cellIdx < U_.size(); ++cellIdx)
    {
        Scalar magSq = U_[cellIdx].magnitudeSquared();
        if (magSq > maxMagSq)
        {
            Scalar mag = std::sqrt(magSq);
            U_[cellIdx] *= maxVelocityMagnitude_ / mag;

            constraintApplied++;
        }
    }

    return constraintApplied;
}

size_t Constraint::applyPressureConstraints() noexcept
{
    if (!enablePressureConstraints_) return 0;

    size_t constraintApplied = 0;

    for (size_t cellIdx = 0; cellIdx < p_.size(); ++cellIdx)
    {
        Scalar clamped = std::clamp(p_[cellIdx], minPressure_, maxPressure_);
        if (clamped != p_[cellIdx])
        {
            p_[cellIdx] = clamped;
            constraintApplied++;
        }
    }

    return constraintApplied;
}
