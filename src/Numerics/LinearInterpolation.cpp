/******************************************************************************
 * @file LinearInterpolation.cpp
 * @brief Implementation of linear interpolation functions
 *****************************************************************************/

#include "LinearInterpolation.hpp"

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

    const Scalar d_P = face.d_Pf_mag();
    const Scalar d_N = face.d_Nf_mag().value();
    const Scalar total = d_P + d_N + vSmallValue;

    const Scalar w_P = d_N / total;
    const Scalar w_N = d_P / total;

    return w_P * field[P] + w_N * field[N];
}

Scalar interpolateToFace
(
    const Face& face,
    const ScalarField& field,
    const BoundaryConditions& bcManager,
    const std::string& fieldName
)
{
    if (face.isBoundary())
    {
        return bcManager.calculateBoundaryFaceValue(face, field, fieldName);
    }

    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();

    const Scalar d_P = face.d_Pf_mag();
    const Scalar d_N = face.d_Nf_mag().value();
    const Scalar total = d_P + d_N + vSmallValue;

    const Scalar w_P = d_N / total;
    const Scalar w_N = d_P / total;

    return w_P * field[P] + w_N * field[N];
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

    const Scalar d_P = face.d_Pf_mag();
    const Scalar d_N = face.d_Nf_mag().value();
    const Scalar total = d_P + d_N + vSmallValue;

    const Scalar w_P = d_N / total;
    const Scalar w_N = d_P / total;

    return w_P * field[P] + w_N * field[N];
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
            face, field, fieldName
        );
    }

    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();

    const Scalar d_P = face.d_Pf_mag();
    const Scalar d_N = face.d_Nf_mag().value();
    const Scalar total = d_P + d_N + vSmallValue;

    const Scalar w_P = d_N / total;
    const Scalar w_N = d_P / total;

    return w_P * field[P] + w_N * field[N];
}