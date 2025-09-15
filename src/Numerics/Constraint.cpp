/******************************************************************************
 * @file Constraint.cpp
 * @brief Implementation of constraint handling for CFD fields
 *****************************************************************************/

#include "Constraint.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

Constraint::Constraint
(
    VectorField& velocityField,
    ScalarField& pressureField
) : U(velocityField),
    p(pressureField),
    enableVelocityConstraints(false),
    enablePressureConstraints(false),
    maxVelocityMagnitude(100.0),         // Default: 1000 m/s
    minPressure(-1e6),                   // Default: -1 MPa
    maxPressure(1e6)                     // Default: 1 MPa
{
}

void Constraint::setVelocityConstraints(Scalar maxVelocity)
{
    this->maxVelocityMagnitude = maxVelocity;
    
    std::cout   << "Velocity constraints set: max velocity = "
                << maxVelocity << " m/s" << std::endl;
}

void Constraint::setPressureConstraints
(
    Scalar minPressure,
    Scalar maxPressure
)
{
    this->minPressure = minPressure;
    this->maxPressure = maxPressure;
    
    std::cout   << "Pressure constraints set: range [" << minPressure
                << ", " << maxPressure << "] Pa" << std::endl;
}

void Constraint::enableConstraints(bool enableVel, bool enablePress)
{
    this->enableVelocityConstraints = enableVel;
    this->enablePressureConstraints = enablePress;
    
    std::cout   << "Field constraints: velocity " 
                << (enableVel ? "enabled" : "disabled")
                << ", pressure " 
                << (enablePress ? "enabled" : "disabled") << std::endl;
}

int Constraint::applyVelocityConstraints()
{
    if (!enableVelocityConstraints) return 0;
    
    int constraintApplications = 0;
    
    for (size_t i = 0; i < U.size(); ++i)
    {
        // Apply velocity magnitude constraint
        Scalar velocityMagnitude = U[i].magnitude();
        if (velocityMagnitude > maxVelocityMagnitude)
        {
            U[i] = U[i] * (maxVelocityMagnitude / velocityMagnitude);
            constraintApplications++;
        }
    }
    
    return constraintApplications;
}

int Constraint::applyPressureConstraints()
{
    if (!enablePressureConstraints) return 0;
    
    int constraintApplications = 0;
    
    for (size_t i = 0; i < p.size(); ++i)
    {
        // Apply pressure bounds
        if (p[i] < minPressure)
        {
            p[i] = minPressure;
            constraintApplications++;
        }
        else if (p[i] > maxPressure)
        {
            p[i] = maxPressure;
            constraintApplications++;
        }
    }
    
    return constraintApplications;
}