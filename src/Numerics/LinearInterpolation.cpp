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
        std::cerr   << "WARNING: interpolateToFace() called on boundary face "
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
        std::cerr   << "WARNING: interpolateToFace() called on boundary face "
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

Vector interpolateToFace
(
    const Face& face,
    const ScalarField& x_component,
    const ScalarField& y_component,
    const ScalarField& z_component,
    const BoundaryConditions& bcManager
)
{
    if (face.isBoundary())
    {
        Scalar vx = 
            bcManager.calculateBoundaryFaceValue(face, x_component, "U_x");
        Scalar vy = 
            bcManager.calculateBoundaryFaceValue(face, y_component, "U_y");
        Scalar vz = 
            bcManager.calculateBoundaryFaceValue(face, z_component, "U_z");
        return Vector(vx, vy, vz);
    }

    const size_t P = face.ownerCell();
    const size_t N = face.neighborCell().value();

    const Scalar d_P = face.d_Pf_mag();
    const Scalar d_N = face.d_Nf_mag().value();
    const Scalar total = d_P + d_N + vSmallValue;

    const Scalar w_P = d_N / total;
    const Scalar w_N = d_P / total;

    Scalar x_f = w_P * x_component[P] + w_N * x_component[N];
    Scalar y_f = w_P * y_component[P] + w_N * y_component[N];
    Scalar z_f = w_P * z_component[P] + w_N * z_component[N];

    return Vector(x_f, y_f, z_f);
}
