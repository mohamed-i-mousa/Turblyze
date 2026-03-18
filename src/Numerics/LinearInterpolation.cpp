/******************************************************************************
 * @file LinearInterpolation.cpp
 * @brief Implementation of linear interpolation functions
 *****************************************************************************/

#include "LinearInterpolation.hpp"
#include "BoundaryConditions.hpp"

#include <iostream>

Scalar interpolateToFace
(
    const Face& face,
    const ScalarField& field
)
{
    if (face.isBoundary())
    {
        std::cerr
            << "WARNING: interpolateToFace() called on boundary face "
            << "without BoundaryConditions. Use the overload with "
            << "bcManager instead. "
            << "Falling back to owner cell value.\n";

        return field[face.ownerCell()];
    }

    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();

    const Scalar dP = face.dPfMag();
    const Scalar dN = face.dNfMag().value();
    const Scalar total = dP + dN + vSmallValue;

    const Scalar wP = dN / total;
    const Scalar wN = dP / total;

    return wP * field[P] + wN * field[N];
}

Vector interpolateToFace
(
    const Face& face,
    const VectorField& field
)
{
    if (face.isBoundary())
    {
        std::cerr
            << "WARNING: interpolateToFace() called on boundary "
            << "face without BoundaryConditions. Use the overload "
            << "with bcManager instead. "
            << "Falling back to owner cell value.\n";

        return field[face.ownerCell()];
    }

    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();

    const Scalar dP = face.dPfMag();
    const Scalar dN = face.dNfMag().value();
    const Scalar total = dP + dN + vSmallValue;

    const Scalar wP = dN / total;
    const Scalar wN = dP / total;

    return wP * field[P] + wN * field[N];
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

    const Scalar dP = face.dPfMag();
    const Scalar dN = face.dNfMag().value();
    const Scalar total = dP + dN + vSmallValue;

    const Scalar wP = dN / total;
    const Scalar wN = dP / total;

    return wP * field[P] + wN * field[N];
}
