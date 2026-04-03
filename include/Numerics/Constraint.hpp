/******************************************************************************
 * @file Constraint.hpp
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

#include "CellData.hpp"


class Constraint
{
public:

    /**
     * @brief Construct constraint manager with field references
     * @param velocityField Reference to velocity field to constrain
     * @param pressureField Reference to pressure field to constrain
     */
    Constraint
    (
        VectorField& velocityField,
        ScalarField& pressureField
    ) noexcept;

    /// Copy and Move Semantics
    Constraint(const Constraint&) = delete;
    Constraint& operator=(const Constraint&) = delete;
    
// Setter methods

    /**
     * @brief Set velocity field constraints
     * @param maxVelocity Maximum allowed velocity magnitude
     */
    void setVelocityConstraints(Scalar maxVelocity) noexcept;

    /**
     * @brief Set pressure field constraints
     * @param minPressure Minimum allowed pressure
     * @param maxPressure Maximum allowed pressure
     */
    void setPressureConstraints(Scalar minPressure, Scalar maxPressure);

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
    [[nodiscard("Constrained cells count needed for diagnostics")]]
    size_t applyVelocityConstraints() noexcept;

    /**
     * @brief Apply pressure field constraints
     * @return Number of cells where constraints were applied
     */
    [[nodiscard("Constrained cells count needed for diagnostics")]]
    size_t applyPressureConstraints() noexcept;

private:

// Private members

    /// Velocity field reference
    VectorField& U_;

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
