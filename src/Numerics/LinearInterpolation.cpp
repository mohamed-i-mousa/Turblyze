/******************************************************************************
 * @file LinearInterpolation.cpp
 * @brief Implementation of linear interpolation functions
 *****************************************************************************/

#include "LinearInterpolation.hpp"
#include "BoundaryConditions.hpp"

#include <cassert>

// Helper function: compute weight for neighbor cell in linear interpolation.
// Returns wN = dP / (dP + dN), so wP = 1 - wN.
// Weight for the closer cell is larger, as expected.
static Scalar faceWeight(const Face& face)
{
    const Scalar dP = face.dPfMag();
    const Scalar dN = face.dNfMag().value();
    return dP / (dP + dN);
}

Scalar interpolateToFace
(
    const Face& face,
    const ScalarField& field
)
{
    assert(!face.isBoundary() && "interpolateToFace (no BC) must not be called on boundary faces");

    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();
    const Scalar wN = faceWeight(face);

    return (S(1.0) - wN) * field[P] + wN * field[N];
}

Vector interpolateToFace
(
    const Face& face,
    const VectorField& field
)
{
    assert(!face.isBoundary() && "interpolateToFace (no BC) must not be called on boundary faces");

    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();
    const Scalar wN = faceWeight(face);

    return (S(1.0) - wN) * field[P] + wN * field[N];
}

Vector interpolateToFace
(
    const Face& face,
    const VectorField& field,
    const BoundaryConditions& bcManager,
    const std::string& fieldName
)
{
    if (face.isBoundary())
    {
        return bcManager.calculateBoundaryVectorFaceValue
        (
            fieldName, field, face
        );
    }

    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();
    const Scalar wN = faceWeight(face);

    return (S(1.0) - wN) * field[P] + wN * field[N];
}
