/******************************************************************************
 * @file LinearInterpolation.hpp
 * @brief Linear interpolation functions for face value reconstruction
 *
 * @details Provides distance-weighted linear interpolation from cell centers 
 * to face centers. Supports both Scalar and Vector fields with optional 
 * BC handling.
 * Uses phif = phiP * wP + phiN * wN with geometric weighting.
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
 * @brief Linear interpolation of vector field with BC handling
 * @param face Face for interpolation
 * @param field Cell-centered vector field
 * @param bcManager Boundary condition manager
 * @param fieldName Name of field for BC lookup (e.g., "U")
 * @return Interpolated vector value with proper BC handling
 */
Vector interpolateToFace
(
    const Face& face,
    const VectorField& field,
    const BoundaryConditions& bcManager,
    const std::string& fieldName
);

#endif // LINEAR_INTERPOLATION_HPP
