#ifndef MASS_FLOW_RATE_H
#define MASS_FLOW_RATE_H

#include <vector>
#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "CellData.h"
#include "FaceData.h"

/**
 * Calculate mass flow rate at each face
 * 
 * For boundary faces: mdot = rho * dot(U_field[P], Sf)
 * For internal faces: mdot = rho * dot([U_field[P] * (1 - w) + U_field[N] * w], Sf)
 * where w is the distance-weighted interpolation factor
 * 
 * @param faces Vector of all faces in the mesh
 * @param cells Vector of all cells in the mesh  
 * @param U_field Velocity field (cell-centered)
 * @param rho Density (constant)
 * @return FaceField containing mass flow rate for each face
 */
FaceFluxField calculateMassFlowRate(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const VectorField& U_field,
    Scalar rho
);

#endif // MASS_FLOW_RATE_H
