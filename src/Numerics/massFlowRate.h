#ifndef MASS_FLOW_RATE_H
#define MASS_FLOW_RATE_H

#include <vector>
#include <map>
#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "CellData.h"
#include "FaceData.h"
#include "BoundaryConditions.h"

/**
 * Calculate mass flow rate at each face
 * 
 * For boundary faces: 
 * - FIXED_VALUE: Use the fixed velocity value from BC
 * - FIXED_GRADIENT: Extrapolate velocity using gradient
 * - ZERO_GRADIENT: Use owner cell velocity (current implementation)
 * - NO_SLIP: Use zero velocity (wall boundary)
 * 
 * For internal faces: mdot = rho * dot([U_field[P] * (1 - w) + U_field[N] * w], Sf)
 * where w is the distance-weighted interpolation factor
 */
FaceFluxField calculateMassFlowRate(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const VectorField& U_field,
    Scalar rho,
    const BoundaryConditions& bcManager,
    const std::map<size_t, const BoundaryPatch*>& faceToPatchMap
);

#endif // MASS_FLOW_RATE_H
