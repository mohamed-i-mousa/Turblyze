/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FixedValue.cpp
 * @brief Fixed-value boundary coefficients
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "FixedValue.h"

// Project headers
#include "BoundaryPatch.h"
#include "Face.h"
#include "Mesh.h"


// ****************************** Public Methods ******************************

void FixedValue::updateCoeffs
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
        const Scalar gDiff = mesh.gDiff(face);

        coeffs_.a[localIdx] = S(0.0);
        coeffs_.b[localIdx] = value_;
        coeffs_.c[localIdx] = -gDiff;
        coeffs_.d[localIdx] = value_ * gDiff;
    }
}

