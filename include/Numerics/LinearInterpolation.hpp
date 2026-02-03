/******************************************************************************
 * @file linearInterpolation.hpp
 * @brief Linear interpolation functions for face value reconstruction
 *
 * Provides distance-weighted linear interpolation from cell centers to face
 * centers. Supports both Scalar and Vector fields with optional BC handling.
 * Uses phi_f = phi_P * w_P + phi_N * w_N with geometric weighting.
 *****************************************************************************/

#ifndef LINEAR_INTERPOLATION_HPP
#define LINEAR_INTERPOLATION_HPP

#include "Face.hpp"
#include "Scalar.hpp"
#include "Vector.hpp"
#include "CellData.hpp"
#include "BoundaryConditions.hpp"

/**
 * @brief Linear interpolation of scalar field to face
 * @param face Face for interpolation
 * @param field Cell-centered scalar field
 * @return Interpolated scalar value at face
 *
 * For boundary faces, returns owner cell value (zero-gradient assumption).
 */
Scalar interpolateToFace
(
    const Face& face,
    const ScalarField& field
);

/**
 * @brief Linear interpolation of scalar field with BC handling
 * @param face Face for interpolation
 * @param field Cell-centered scalar field
 * @param bcManager Boundary condition manager
 * @param fieldName Name of field for BC lookup (e.g., "p", "k")
 * @return Interpolated scalar value with proper BC handling
 */
Scalar interpolateToFace
(
    const Face& face,
    const ScalarField& field,
    const BoundaryConditions& bcManager,
    const std::string& fieldName
);

/**
 * @brief Linear interpolation of vector field to face
 * @param face Face for interpolation
 * @param field Cell-centered vector field
 * @return Interpolated vector value at face
 *
 * For boundary faces, returns owner cell value (zero-gradient assumption).
 */
Vector interpolateToFace
(
    const Face& face,
    const VectorField& field
);

/**
 * @brief Linear interpolation of vector from component fields with BC handling
 * @param face Face for interpolation
 * @param x_component X-component of vector field
 * @param y_component Y-component of vector field
 * @param z_component Z-component of vector field
 * @param bcManager Boundary condition manager
 * @return Interpolated vector value with proper BC handling
 *
 * Reconstructs a vector field from its scalar components with BC evaluation.
 * Each component uses scalar BC method at boundaries.
 */
Vector interpolateToFace
(
    const Face& face,
    const ScalarField& x_component,
    const ScalarField& y_component,
    const ScalarField& z_component,
    const BoundaryConditions& bcManager
);

#endif // LINEAR_INTERPOLATION_HPP
