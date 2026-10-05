/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Symmetry.cpp
 * @brief Symmetry-plane boundary coefficients and registration
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Symmetry.h"

// Standard library headers
#include <algorithm>

// Project headers
#include "BoundaryPatch.h"
#include "Face.h"
#include "Mesh.h"
#include "Vector.h"

// ***************************** Internal Helpers *****************************

namespace
{

// Entry of a vector along the velocity component
Scalar along(const Vector& v, Field component) noexcept
{
    switch (component)
    {
        case Field::Ux:
            return v.x();
        case Field::Uy:
            return v.y();
        case Field::Uz:
            return v.z();
        default:
            return S(0.0);
    }
}

// Sum of n_j U_j over the two components other than the given one
Scalar crossNormalVelocity
(
    const Vector& n,
    const Vector& U,
    Field component
) noexcept
{
    switch (component)
    {
        case Field::Ux:
            return n.y() * U.y() + n.z() * U.z();
        case Field::Uy:
            return n.x() * U.x() + n.z() * U.z();
        case Field::Uz:
            return n.x() * U.x() + n.y() * U.y();
        default:
            return S(0.0);
    }
}

} // namespace

// ****************************** Public Methods ******************************

void Symmetry::updateCoeffs
(
    const Mesh& mesh,
    const BoundaryPatch& patch
)
{
    coeffs_ = BoundaryCoeffs(patch.numFaces());

    const Count numFaces = patch.numFaces();

    for (Index localIdx = 0; localIdx < numFaces; ++localIdx)
    {
        const Index faceIdx = patch.firstFaceIdx() + localIdx;
        const Face& face = mesh.faces()[faceIdx];
        const Scalar ni = along(face.normal(), component_);
        const Scalar gDiff = mesh.gDiff(face);

        coeffs_.a[localIdx] = S(1.0) - ni * ni;
        coeffs_.b[localIdx] = S(0.0);
        coeffs_.c[localIdx] = -ni * ni * gDiff;
        coeffs_.d[localIdx] = S(0.0);
    }
}


void Symmetry::refreshCoeffs
(
    const Mesh& mesh,
    const BoundaryPatch& patch,
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz
)
{
    const Count numFaces = patch.numFaces();

    for (Index localIdx = 0; localIdx < numFaces; ++localIdx)
    {
        const Index faceIdx = patch.firstFaceIdx() + localIdx;
        const Face& face = mesh.faces()[faceIdx];
        const Vector n = face.normal();
        const Index owner = face.ownerCell();

        const Vector Uowner(Ux[owner], Uy[owner], Uz[owner]);
        const Scalar ni = along(n, component_);
        const Scalar unCross = crossNormalVelocity(n, Uowner, component_);
        const Scalar gDiff = mesh.gDiff(face);

        coeffs_.b[localIdx] = -ni * unCross;
        coeffs_.d[localIdx] = -ni * unCross * gDiff;
    }
}
