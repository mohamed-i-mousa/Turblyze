#ifndef LINEAR_INTERPOLATION_H
#define LINEAR_INTERPOLATION_H

#include "Face.h"
#include "Scalar.h"
#include "Vector.h"
#include "CellData.h"
#include "BoundaryConditions.h"

/**
 * @brief Compute linear interpolation weights for a face
 * @param face Face for interpolation
 * @param w_P Weight for owner cell
 * @param w_N Weight for neighbor cell
 */
void computeLinearWeights
(
    const Face& face,
    Scalar& w_P,
    Scalar& w_N
);

/**
 * @brief Linear interpolation of vector field to face
 * @param face Face for interpolation
 * @param cellField Cell-centered vector field
 * @return Interpolated vector value at face
 */
Vector VectorLinearInterpolation
(
    const Face& face,
    const VectorField& cellField
);

/**
 * @brief Linear interpolation of vector field to face with boundary condition support
 * @param face Face for interpolation
 * @param cellField Cell-centered vector field
 * @param bcManager Boundary condition manager
 * @param fieldName Field name for boundary condition lookup
 * @return Interpolated vector value at face with proper BC handling
 */
Vector VectorLinearInterpolation
(
    const Face& face,
    const VectorField& cellField,
    const BoundaryConditions& bcManager,
    const std::string& fieldName
);

/**
 * @brief Linear interpolation of scalar field to face
 * @param face Face for interpolation
 * @param cellField Cell-centered scalar field
 * @return Interpolated scalar value at face
 */
Scalar linearInterpolation
(
    const Face& face,
    const ScalarField& cellField
);

#endif