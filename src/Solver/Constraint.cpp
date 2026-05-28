/******************************************************************************
 * @file Constraint.cpp
 * @brief Implementation of constraint handling for CFD fields
 *****************************************************************************/

#include "Constraint.h"

#include <cmath>

Constraint::Constraint
(
    ScalarField& Ux,
    ScalarField& Uy,
    ScalarField& Uz,
    ScalarField& pressureField,
    const bool velocityEnabled,
    const bool pressureEnabled,
    const Scalar maxVelocityMagnitude,
    const Scalar minPressure,
    const Scalar maxPressure
) noexcept
:
    Ux_(Ux),
    Uy_(Uy),
    Uz_(Uz),
    p_(pressureField),
    enableVelocityConstraints_(velocityEnabled),
    enablePressureConstraints_(pressureEnabled),
    maxVelocityMagnitude_(maxVelocityMagnitude),
    minPressure_(minPressure),
    maxPressure_(maxPressure)
{}

size_t Constraint::applyVelocityConstraints() noexcept
{
    if (!enableVelocityConstraints_) return 0;

    size_t constraintApplied = 0;
    const Scalar maxMagSq = maxVelocityMagnitude_ * maxVelocityMagnitude_;

    for (size_t cellIdx = 0; cellIdx < Ux_.size(); ++cellIdx)
    {
        const Scalar magSq =
            Ux_[cellIdx] * Ux_[cellIdx]
          + Uy_[cellIdx] * Uy_[cellIdx]
          + Uz_[cellIdx] * Uz_[cellIdx];

        if (magSq > maxMagSq)
        {
            const Scalar scale =
                maxVelocityMagnitude_ / std::sqrt(magSq);

            Ux_[cellIdx] *= scale;
            Uy_[cellIdx] *= scale;
            Uz_[cellIdx] *= scale;

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
        const Scalar clamped =
            std::clamp(p_[cellIdx], minPressure_, maxPressure_);
            
        if (clamped != p_[cellIdx])
        {
            p_[cellIdx] = clamped;
            constraintApplied++;
        }
    }

    return constraintApplied;
}
