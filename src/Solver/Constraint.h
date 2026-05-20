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

#include "CellData.h"


class Constraint
{
public:

    /**
     * @brief Construct constraint manager with field references
     * @param Ux Reference to x-velocity component to constrain
     * @param Uy Reference to y-velocity component to constrain
     * @param Uz Reference to z-velocity component to constrain
     * @param pressureField Reference to pressure field to constrain
     */
    Constraint
    (
        ScalarField& Ux,
        ScalarField& Uy,
        ScalarField& Uz,
        ScalarField& pressureField
    ) noexcept;

    /// Copy constructor and assignment - Not copyable (reference member)
    Constraint(const Constraint&) = delete;
    Constraint& operator=(const Constraint&) = delete;

    /// Move constructor and assignment - Not movable (reference member)
    Constraint(Constraint&&) = delete;
    Constraint& operator=(Constraint&&) = delete;

    /// Destructor
    ~Constraint() noexcept = default;

// Setter methods

    /**
     * @brief Set velocity field constraints
     * @param maxVelocity Maximum allowed velocity magnitude
     */
    void setVelocityConstraints(Scalar maxVelocity) noexcept
    {
        maxVelocityMagnitude_ = std::abs(maxVelocity);
    }

    /**
     * @brief Set pressure field constraints
     * @param minPressure Minimum allowed pressure
     * @param maxPressure Maximum allowed pressure
     */
    void setPressureConstraints
    (
        Scalar minPressure,
        Scalar maxPressure
    ) noexcept;

    /**
     * @brief Enable or disable field constraints
     * @param enableVel Enable velocity constraints
     * @param enablePress Enable pressure constraints
     */
    void enableConstraints
    (
        bool enableVel,
        bool enablePress
    ) noexcept;

    /**
     * @brief Apply velocity field constraints
     * @return Number of cells where constraints were applied
     */
    [[nodiscard]] size_t applyVelocityConstraints() noexcept;

    /**
     * @brief Apply pressure field constraints
     * @return Number of cells where constraints were applied
     */
    [[nodiscard]] size_t applyPressureConstraints() noexcept;

private:

// Private members

    /// x-velocity component reference
    ScalarField& Ux_;

    /// y-velocity component reference
    ScalarField& Uy_;

    /// z-velocity component reference
    ScalarField& Uz_;

    /// Pressure field reference
    ScalarField& p_;

    /// Enable velocity field constraints
    bool enableVelocityConstraints_ = false;

    /// Enable pressure field constraints
    bool enablePressureConstraints_ = false;

    /// Maximum allowed velocity magnitude
    Scalar maxVelocityMagnitude_ = S(100);

    /// Minimum allowed pressure
    Scalar minPressure_ = S(-1e6);

    /// Maximum allowed pressure
    Scalar maxPressure_ = S(1e6);
};
