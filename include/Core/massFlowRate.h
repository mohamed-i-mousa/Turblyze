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
 * @brief Calculate mass flow rate at each face
 * 
 * @param faces Vector of mesh faces
 * @param cells Vector of mesh cells  
 * @param U_field Velocity field at cell centers
 * @param bcManager Boundary conditions manager
 * @param faceToPatchMap Mapping from face ID to boundary patch
 * @return FaceFluxField Mass flow rate at each face
 * 
 * For internal faces: 
 *      mdot = rho * dot([U_field[P] * (1 - w) + U_field[N] * w], Sf)
 *      where w is the distance-weighted interpolation factor.
 * 
 * For boundary faces, uses boundary condition handling from bcManager.
 * 
 */
FaceFluxField calculateMassFlowRate
(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const VectorField& U_field,
    const BoundaryConditions& bcManager,
    const std::map<size_t, const BoundaryPatch*>& faceToPatchMap
);

#endif