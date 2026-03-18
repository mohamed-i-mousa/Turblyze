/******************************************************************************
 * @file Constraint.cpp
 * @brief Implementation of constraint handling for CFD fields
 *****************************************************************************/

#include "Constraint.hpp"

#include <cmath>

Constraint::Constraint
(
    VectorField& velocityField,
    ScalarField& pressureField
) noexcept
  : U_(velocityField),
    p_(pressureField)
{
}

void Constraint::setVelocityConstraints(Scalar maxVelocity) noexcept
{
    maxVelocityMagnitude_ = maxVelocity;
}

void Constraint::setPressureConstraints
(
    Scalar minPressure,
    Scalar maxPressure
) noexcept
{
    minPressure_ = minPressure;
    maxPressure_ = maxPressure;
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

    size_t constraintApplications = 0;
    Scalar maxMagSq = maxVelocityMagnitude_ * maxVelocityMagnitude_;

    for (size_t cellIdx = 0; cellIdx < U_.size(); ++cellIdx)
    {
        Scalar magSq = U_[cellIdx].magnitudeSquared();
        if (magSq > maxMagSq)
        {
            Scalar mag = std::sqrt(magSq);
            U_[cellIdx] *= maxVelocityMagnitude_ / mag;

            constraintApplications++;
        }
    }

    return constraintApplications;
}

size_t Constraint::applyPressureConstraints() noexcept
{
    if (!enablePressureConstraints_) return 0;

    size_t constraintApplications = 0;

    for (size_t cellIdx = 0; cellIdx < p_.size(); ++cellIdx)
    {
        // Apply pressure bounds
        if (p_[cellIdx] < minPressure_)
        {
            p_[cellIdx] = minPressure_;
            constraintApplications++;
        }
        else if (p_[cellIdx] > maxPressure_)
        {
            p_[cellIdx] = maxPressure_;
            constraintApplications++;
        }
    }

    return constraintApplications;
}
