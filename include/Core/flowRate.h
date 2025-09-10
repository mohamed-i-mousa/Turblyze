/******************************************************************************
 * @file flowRate.h
 * @brief Mass flow rate calculation for finite volume faces
 * 
 * Calculates volumetric flow rates at face centers using linear interpolation
 * of velocity fields with boundary condition integration. For internal faces,
 * uses distance-weighted interpolation: mdot = dot(U_interp, Sf).
 * For boundary faces, applies proper boundary condition handling.
 *****************************************************************************/

#ifndef FLOW_RATE_H
#define FLOW_RATE_H

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
 * @brief Calculate volumetric flow rate at each face
 * @param faces Vector of mesh faces
 * @param U_field Velocity field at cell centers
 * @param bcManager Boundary conditions manager
 * @param faceToPatchMap Mapping from face ID to boundary patch
 * @return FaceFluxField volumetric flow rate at each face
 * 
 * For internal faces: mDot = Dot([U_field[P] * (1-w) + U_field[N] * w], Sf)
 * where w is the distance-weighted interpolation factor.
 * For boundary faces, uses boundary condition handling from bcManager.
 * 
 * @note Incompressible equations are divided by the density.
 *       Then, mDot is the volumetric flow rate.
 */
FaceFluxField calculateFlowRate
(
    const std::vector<Face>& faces,
    const VectorField& U_field,
    const BoundaryConditions& bcManager,
    const std::map<size_t, const BoundaryPatch*>& faceToPatchMap
);

#endif