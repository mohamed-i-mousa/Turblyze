/******************************************************************************
 * @file Constraint.h
 * @brief Field constraint system for CFD solution stability
 *
 * @details This header defines the Constraint class to apply constraints on
 * velocity and pressure fields to prevent overshooting and divergence.
 *
 * @class Constraint
 *
 * The Constraint class manages:
 * - Velocity magnitude clipping (preventing unrealistically high speeds)
 * - Pressure bounds (preventing excessive pressures)
 * - Reporting of constraint activity (number of cells clipped)
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "CellData.h"

// ***************************** class Constraint *****************************

class Constraint
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    Constraint
    (
        ScalarField& Ux,
        ScalarField& Uy,
        ScalarField& Uz,
        ScalarField& pressureField,
        const bool velocityEnabled,
        const bool pressureEnabled,
        const Scalar maxVelocityMagnitude,
        const Scalar minPressure,
        const Scalar maxPressure,
        const bool debug
    ) noexcept;

    /// Copy constructor and assignment - Not copyable (reference member)
    Constraint(const Constraint&) = delete;
    Constraint& operator=(const Constraint&) = delete;

    /// Move constructor and assignment - Not movable (reference member)
    Constraint(Constraint&&) = delete;
    Constraint& operator=(Constraint&&) = delete;

    /// Destructor
    ~Constraint() noexcept = default;

// ****************************** Public Methods ******************************

    /// Apply velocity field constraints
    void applyVelocityConstraints() noexcept;

    /// Apply pressure field constraints
    void applyPressureConstraints() noexcept;

// ****************************** Private Members *****************************

private:

    /// x-velocity component reference
    ScalarField& Ux_;

    /// y-velocity component reference
    ScalarField& Uy_;

    /// z-velocity component reference
    ScalarField& Uz_;

    /// Pressure field reference
    ScalarField& p_;

    /// Enable velocity field constraints
    bool enableVelocityConstraints_;

    /// Enable pressure field constraints
    bool enablePressureConstraints_;

    /// Maximum allowed velocity magnitude
    Scalar maxVelocityMagnitude_;

    /// Minimum allowed pressure
    Scalar minPressure_;

    /// Maximum allowed pressure
    Scalar maxPressure_;

    /// Report when constraints clip cells
    bool debug_;
};
