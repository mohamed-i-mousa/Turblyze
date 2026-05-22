/******************************************************************************
 * @file LinearInterpolation.h
 * @brief Linear interpolation of a cell-centered field to internal faces
 *
 * @details Distance-weighted linear interpolation of CellData<T> from owner
 * and neighbour cell centres to the shared face. For Boundary-face values,
 * use BoundaryConditions::boundaryFaceValue before calling this function.
 *****************************************************************************/

#pragma once

#include "Face.h"
#include "CellData.h"
#include "ErrorHandler.h"

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
    const Face& face,
    const CellData<T>& field
)
{
    if (face.isBoundary())
    {
        FatalError
        (
            "interpolateToFace must not be called on boundary "
            "faces; resolve BC values at the call site"
        );
    }

    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();
    const Scalar wN = faceWeight(face);

    return (S(1.0) - wN) * field[P] + wN * field[N];
}
