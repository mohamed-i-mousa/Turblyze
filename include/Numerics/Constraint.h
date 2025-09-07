/******************************************************************************
 * @file Constraint.h
 * @brief Field constraint system for CFD solution stability
 * 
 * This class provides mechanisms to apply constraints on velocity and pressure
 * fields to prevent overshooting and divergence during CFD iterations.
 * The constraints include magnitude limits and value bounds
 * to maintain numerical stability.
 * 
 * @author Mohamed Mousa
 * @date 2025
 *****************************************************************************/

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "Scalar.h"
#include "Vector.h"
#include "CellData.h"

/**
 * @class Constraint
 * @brief Manages field constraints for CFD solution stability
 * 
 * This class applies various constraints on velocity and pressure fields
 * to prevent numerical instabilities:
 * - Velocity magnitude limits
 * - Pressure bound limits
 */
class Constraint
{
public:
    /**
     * @brief Constructor
     * @param velocityField Reference to velocity field
     * @param pressureField Reference to pressure field
     */
    Constraint(
        VectorField& velocityField,
        ScalarField& pressureField
    );
    
    /**
     * @brief Set velocity field constraints
     * @param maxVelocity Maximum allowed velocity magnitude
     */
    void setVelocityConstraints(Scalar maxVelocity);
    
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
    void enableConstraints(bool enableVel = true, bool enablePress = true);
    
    /**
     * @brief Apply velocity field constraints
     * @return Number of cells where constraints were applied
     */
    int applyVelocityConstraints();
    
    /**
     * @brief Apply pressure field constraints  
     * @return Number of cells where constraints were applied
     */
    int applyPressureConstraints();
    
    /**
     * @brief Check if velocity constraints are enabled
     * @return True if velocity constraints are enabled
     */
    bool areVelocityConstraintsEnabled() const { return enableVelocityConstraints; }
    
    /**
     * @brief Check if pressure constraints are enabled
     * @return True if pressure constraints are enabled
     */
    bool arePressureConstraintsEnabled() const { return enablePressureConstraints; }

private:
    /// Field references
    VectorField& U;                        ///< Velocity field reference
    ScalarField& p;                        ///< Pressure field reference
    
    /// Constraint enable flags
    bool enableVelocityConstraints;        ///< Enable velocity field constraints
    bool enablePressureConstraints;        ///< Enable pressure field constraints
    
    /// Velocity constraint parameters
    Scalar maxVelocityMagnitude;           ///< Maximum allowed velocity magnitude
    
    /// Pressure constraint parameters
    Scalar minPressure;                    ///< Minimum allowed pressure
    Scalar maxPressure;                    ///< Maximum allowed pressure
};

#endif // CONSTRAINT_H