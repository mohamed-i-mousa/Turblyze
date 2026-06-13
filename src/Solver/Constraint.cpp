/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Constraint.cpp
 * @brief Implementation of constraint handling for CFD fields
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Constraint.h"

// Standard library headers
#include <cmath>
#include <iostream>

// ************************* Special Member Functions *************************

Constraint::Constraint
(
    ScalarField& Ux,
    ScalarField& Uy,
    ScalarField& Uz,
    ScalarField& pressureField,
    bool velocityEnabled,
    bool pressureEnabled,
    Scalar maxVelocityMagnitude,
    Scalar minPressure,
    Scalar maxPressure,
    bool debug
) noexcept
:
    Ux_{Ux},
    Uy_{Uy},
    Uz_{Uz},
    p_{pressureField},
    enableVelocityConstraints_{velocityEnabled},
    enablePressureConstraints_{pressureEnabled},
    maxVelocityMagnitude_{maxVelocityMagnitude},
    minPressure_{minPressure},
    maxPressure_{maxPressure},
    debug_{debug}
{}

// ****************************** Public Methods ******************************

void Constraint::applyVelocityConstraints() noexcept
{
    if (!enableVelocityConstraints_) return;

    Count constraintApplied = 0;
    const Scalar maxMagSq = maxVelocityMagnitude_ * maxVelocityMagnitude_;

    #pragma omp parallel for schedule(static) reduction(+:constraintApplied)
    for (Index cellIdx = 0; cellIdx < Ux_.size(); ++cellIdx)
    {
        const Scalar magSq =
            Ux_[cellIdx] * Ux_[cellIdx]
          + Uy_[cellIdx] * Uy_[cellIdx]
          + Uz_[cellIdx] * Uz_[cellIdx];

        if (magSq > maxMagSq)
        {
            const Scalar scale = maxVelocityMagnitude_ / std::sqrt(magSq);

            Ux_[cellIdx] *= scale;
            Uy_[cellIdx] *= scale;
            Uz_[cellIdx] *= scale;

            constraintApplied++;
        }
    }

    if (debug_ && constraintApplied > 0)
    {
        std::cout
            << "  Applied velocity constraints to "
            << constraintApplied << " cells" << '\n';
    }
}

void Constraint::applyPressureConstraints() noexcept
{
    if (!enablePressureConstraints_) return;

    Count constraintApplied = 0;

    #pragma omp parallel for schedule(static) reduction(+:constraintApplied)
    for (Index cellIdx = 0; cellIdx < p_.size(); ++cellIdx)
    {
        const Scalar clamped =
            std::clamp(p_[cellIdx], minPressure_, maxPressure_);

        if (clamped != p_[cellIdx])
        {
            p_[cellIdx] = clamped;
            constraintApplied++;
        }
    }

    if (debug_ && constraintApplied > 0)
    {
        std::cout
            << "  Applied pressure constraints to "
            << constraintApplied << " cells" << '\n';
    }
}
