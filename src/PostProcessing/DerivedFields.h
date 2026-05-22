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

/// Compute velocity magnitude field from velocity components
[[nodiscard]] ScalarField velocityMagnitude
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz
);

/// Compute vorticity magnitude field from vorticity vector field
[[nodiscard]] ScalarField vorticityMagnitude
(
    const VectorField& vorticity
);

/// Compute Q-criterion for vortex identification
/// Q = 0.5 * (||Omega||^2 - ||S||^2)
[[nodiscard]] ScalarField QCriterion
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
);

/// Compute strain rate magnitude = sqrt(2 * S_ij * S_ij)
[[nodiscard]] ScalarField strainRateMagnitude
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
);

} // namespace VTK
