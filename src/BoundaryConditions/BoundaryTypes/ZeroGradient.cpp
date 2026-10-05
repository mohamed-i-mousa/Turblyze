/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ZeroGradient.cpp
 * @brief Zero-gradient boundary coefficients
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "ZeroGradient.h"

// Project headers
#include "BoundaryPatch.h"


// ****************************** Public Methods ******************************

void ZeroGradient::updateCoeffs
(
    const Mesh& /* mesh */,
    const BoundaryPatch& patch
)
{
    // Constructor's default: a = 1, b = c = d = 0
    coeffs_ = BoundaryCoeffs(patch.numFaces());
}

