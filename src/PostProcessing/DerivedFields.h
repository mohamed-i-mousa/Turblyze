/******************************************************************************
 * @file DerivedFields.h
 * @brief Cell-centered derived scalar fields for post-processing
 *
 * @details Pure field-math utilities — no VTK dependencies. Produces scalar
 * fields (magnitudes, Q-criterion, strain rate) from velocity or gradient
 * vector fields. Intended for use ahead of VTK export so consumers can
 * write the derived quantities into `.vtu` files.
 *****************************************************************************/

#pragma once

#include "CellData.h"


namespace VTK
{

/**
 * @brief Compute velocity magnitude field from velocity vector field
 * @param velocity Input velocity vector field
 * @return Scalar field containing velocity magnitude at each cell
 */
[[nodiscard]] ScalarField velocityMagnitude
(
    const VectorField& velocity
);

/**
 * @brief Compute vorticity magnitude field from vorticity vector field
 * @param vorticity Input vorticity vector field
 * @return Scalar field containing vorticity magnitude at each cell
 */
[[nodiscard]] ScalarField vorticityMagnitude
(
    const VectorField& vorticity
);

/**
 * @brief Compute Q-criterion for vortex identification
 *
 * @details Q-criterion identifies vortex cores as regions where Q > 0,
 * where Q = 0.5 * (||Omega||^2 - ||S||^2), with Omega being the rotation
 * rate tensor and S the strain rate tensor.
 *
 * @param gradUx Gradient of x-velocity component
 * @param gradUy Gradient of y-velocity component
 * @param gradUz Gradient of z-velocity component
 * @return Scalar field containing Q-criterion at each cell
 */
[[nodiscard]] ScalarField QCriterion
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
);

/**
 * @brief Compute strain rate magnitude field
 *
 * @details Strain rate magnitude = sqrt(2 * S_ij * S_ij) where S_ij is the
 * symmetric strain rate tensor.
 *
 * @param gradUx Gradient of x-velocity component
 * @param gradUy Gradient of y-velocity component
 * @param gradUz Gradient of z-velocity component
 * @return Scalar field containing strain rate magnitude at each cell
 */
[[nodiscard]] ScalarField strainRateMagnitude
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
);

} // namespace VTK
