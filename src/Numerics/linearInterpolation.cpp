/******************************************************************************
 * @file linearInterpolation.cpp
 * @brief Implementation of linear interpolation functions for CFD fields
 *****************************************************************************/

#include "linearInterpolation.h"

#include <iostream>
#include <cmath>
#include <algorithm>

// ----- Linear interpolation helpers ----- //
void computeLinearWeights
(
    const Face& face,
    Scalar& w_P,
    Scalar& w_N
)
{
    if (face.isBoundary())
    {
        w_P = S(1.0);
        w_N = 0.0;
        return;
    }

    const Scalar d_P = face.d_Pf_mag();
    const Scalar d_N = face.d_Nf_mag().value();
    const Scalar total = d_P + d_N;

    w_P = d_N / (total + vSmallValue);
    w_N = d_P / (total + vSmallValue);

}

Vector VectorLinearInterpolation
(
    const Face& face,
    const VectorField& cellField
)
{
    if (face.isBoundary())
    {
        return cellField[face.ownerCell()];
    }

    Scalar w_P = 0.0, w_N = 0.0;
    computeLinearWeights(face, w_P, w_N);
    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();
    return w_P * cellField[P] + w_N * cellField[N];
}

Vector VectorLinearInterpolation
(
    const Face& face,
    const VectorField& cellField,
    const BoundaryConditions& bcManager,
    const std::string& fieldName
)
{
    if (face.isBoundary())
    {
        return bcManager.calculateBoundaryFaceVectorValue(face, cellField, fieldName);
    }

    Scalar w_P = 0.0, w_N = 0.0;
    computeLinearWeights(face, w_P, w_N);
    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();
    return w_P * cellField[P] + w_N * cellField[N];
}

Scalar linearInterpolation
(
    const Face& face,
    const ScalarField& cellField
)
{
    if (face.isBoundary())
    {
        return cellField[face.ownerCell()];
    }

    Scalar w_P = 0.0, w_N = 0.0;
    computeLinearWeights(face, w_P, w_N);
    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();
    return w_P * cellField[P] + w_N * cellField[N];
}