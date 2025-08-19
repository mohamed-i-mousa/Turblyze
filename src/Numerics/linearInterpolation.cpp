#include "linearInterpolation.h"

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

    const Scalar d_P = face.d_Pf_mag;
    const Scalar d_N = face.d_Nf_mag.value();
    const Scalar total = d_P + d_N + 1e-20;
    w_P = d_N / total;
    w_N = d_P / total;
}

Vector VectorLinearInterpolation
(
    const Face& face,
    const VectorField& cellField
)
{
    if (face.isBoundary())
    {
        return cellField[face.ownerCell];
    }

    Scalar w_P = 0.0, w_N = 0.0;
    computeLinearWeights(face, w_P, w_N);
    const size_t P = face.ownerCell;
    const size_t N = face.neighbourCell.value();
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
        return cellField[face.ownerCell];
    }

    Scalar w_P = 0.0, w_N = 0.0;
    computeLinearWeights(face, w_P, w_N);
    const size_t P = face.ownerCell;
    const size_t N = face.neighbourCell.value();
    return w_P * cellField[P] + w_N * cellField[N];
}