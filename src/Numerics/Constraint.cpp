/******************************************************************************
 * @file Constraint.cpp
 * @brief Implementation of constraint handling for CFD fields
 *****************************************************************************/

#include "Constraint.hpp"

Constraint::Constraint
(
    VectorField& velocityField,
    ScalarField& pressureField
) : U_(velocityField),
    p_(pressureField),
    enableVelocityConstraints_(false),
    enablePressureConstraints_(false),
    maxVelocityMagnitude_(100.0),
    minPressure_(-1e6),
    maxPressure_(1e6)          
{
}

void Constraint::setVelocityConstraints(Scalar maxVelocity)
{
    maxVelocityMagnitude_ = maxVelocity;
}

void Constraint::setPressureConstraints
(
    Scalar minPressure,
    Scalar maxPressure
)
{
    minPressure_ = minPressure;
    maxPressure_ = maxPressure;
}

void Constraint::enableConstraints(bool enableVel, bool enablePress)
{
    enableVelocityConstraints_ = enableVel;
    enablePressureConstraints_ = enablePress;
}

size_t Constraint::applyVelocityConstraints()
{
    if (!enableVelocityConstraints_) return 0;

    size_t constraintApplications = 0;

    for (size_t CellIdx = 0; CellIdx < U_.size(); ++CellIdx)
    {
        Scalar velocityMagnitude = U_[CellIdx].magnitude();
        if (velocityMagnitude > maxVelocityMagnitude_)
        {
            U_[CellIdx] = 
                U_[CellIdx] * (maxVelocityMagnitude_ / velocityMagnitude);

            constraintApplications++;
        }
    }

    return constraintApplications;
}

size_t Constraint::applyPressureConstraints()
{
    if (!enablePressureConstraints_) return 0;

    size_t constraintApplications = 0;

    for (size_t CellIdx = 0; CellIdx < p_.size(); ++CellIdx)
    {
        // Apply pressure bounds
        if (p_[CellIdx] < minPressure_)
        {
            p_[CellIdx] = minPressure_;
            constraintApplications++;
        }
        else if (p_[CellIdx] > maxPressure_)
        {
            p_[CellIdx] = maxPressure_;
            constraintApplications++;
        }
    }

    return constraintApplications;
}