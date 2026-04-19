/******************************************************************************
 * @file LinearInterpolation.hpp
 * @brief Linear interpolation functions for face value reconstruction
 *
 * @details Provides distance-weighted linear interpolation from cell centers
 * to face centers. Supports both Scalar and Vector fields with optional
 * BC handling.
 * Uses phif = phiP * wP + phiN * wN with geometric weighting.
 *****************************************************************************/

#pragma once

#include "Face.hpp"
#include "CellData.hpp"

/// Forward declaration
class BoundaryConditions;

/**
 * @brief Linear interpolation of scalar field to face
 * @param face Face for interpolation
 * @param field Cell-centered scalar field
 * @return Interpolated scalar value at face
 */
[[nodiscard]] Scalar interpolateToFace
(
    const Face& face,
    const ScalarField& field
);

/**
 * @brief Linear interpolation of vector field to face
 * @param face Face for interpolation
 * @param field Cell-centered vector field
 * @return Interpolated vector value at face
 *
 * For boundary faces, returns owner cell value (zero-gradient assumption).
 */
[[nodiscard]] Vector interpolateToFace
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
[[nodiscard]] Vector interpolateToFace
(
    const Face& face,
    const VectorField& field,
    const BoundaryConditions& bcManager,
    const std::string& fieldName
);

/**
 * @brief Linear interpolation of tensor field to face
 * @param face Face for interpolation
 * @param field Cell-centered tensor field
 * @return Interpolated tensor value at face
 *
 * For boundary faces, callers must handle BCs explicitly; this
 * overload terminates the program if called on a boundary face.
 */
[[nodiscard]] Tensor interpolateToFace
(
    const Face& face,
    const TensorField& field
);
