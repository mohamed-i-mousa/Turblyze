/******************************************************************************
 * @file LinearInterpolation.hpp
 * @brief Linear interpolation of a cell-centered field to internal faces
 *
 * @details Distance-weighted linear interpolation of CellData<T> from owner
 * and neighbour cell centres to the shared face. Boundary-face values are
 * the caller's responsibility — use BoundaryConditions::boundaryVectorFaceValue
 * (or the analogous helper for scalar/tensor fields) before invoking this
 * function.
 *****************************************************************************/

#pragma once

#include "Face.hpp"
#include "CellData.hpp"
#include "ErrorHandler.hpp"

/// Distance weight for the neighbour cell: wN = dP / (dP + dN).
[[nodiscard]] inline Scalar faceWeight(const Face& face)
{
    const Scalar dP = face.dPfMag();
    const Scalar dN = face.dNfMag().value();
    return dP / (dP + dN);
}

/**
 * @brief Linear interpolation of a cell-centered field to an internal face
 * @param face Internal face (FatalError if boundary)
 * @param field Cell-centered field of type Scalar, Vector, or Tensor
 * @return Face-centered value
 * @note Callers must resolve boundary-face values themselves.
 */
template<typename T>
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
