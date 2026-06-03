/******************************************************************************
 * @file LinearInterpolation.h
 * @brief Linear interpolation of a cell-centered field to internal faces
 *
 * @details Distance-weighted linear interpolation of CellData<T> from owner
 * and neighbour cell centres to the shared face. For Boundary-face values,
 * use BoundaryConditions::boundaryFaceValue before calling this function.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Face.h"
#include "CellData.h"
#include "ErrorHandler.h"

// ************************** Interpolation Functions *************************

/// Distance weight for the neighbour cell: wN = dP / (dP + dN).
[[nodiscard]] inline Scalar faceWeight(const Face& face)
{
    const Scalar dP = face.dPfMag();
    const Scalar dN = face.dNfMag().value();
    return dP / (dP + dN);
}

/// Linear interpolation of a cell-centered field to an internal face
template<CellFieldType T>
[[nodiscard]] T interpolateToFace
(
    const Face& targetFace,
    const CellData<T>& field
)
{
    if (targetFace.isBoundary())
    {
        FatalError
        (
            "interpolateToFace must not be called on boundary "
            "faces; resolve BC values at the call site"
        );
    }

    const Index P = targetFace.ownerCell();
    const Index N = targetFace.neighborCell().value();
    const Scalar wN = faceWeight(targetFace);

    return (S(1.0) - wN) * field[P] + wN * field[N];
}
