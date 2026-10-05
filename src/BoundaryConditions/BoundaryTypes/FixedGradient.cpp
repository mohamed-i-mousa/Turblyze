/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FixedGradient.cpp
 * @brief Fixed-gradient boundary coefficients and diagnostics
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "FixedGradient.h"

// Project headers
#include "BoundaryPatch.h"
#include "Face.h"
#include "Mesh.h"


// ****************************** Public Methods ******************************

void FixedGradient::updateCoeffs
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
        const Scalar normalDistance = dot(mesh.dPf(face), face.normal());

        coeffs_.a[localIdx] = S(1.0);
        coeffs_.b[localIdx] = gradient_ * normalDistance;
        coeffs_.c[localIdx] = S(0.0);
        coeffs_.d[localIdx] = gradient_;
    }
}

