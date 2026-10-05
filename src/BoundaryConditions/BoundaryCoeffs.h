/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryCoeffs.h
 * @brief Linearized boundary condition coefficient containers
 *
 * @details Boundary conditions are linearized into 4 coefficients
 * a and b are the value coefficients that construct the face value as:
 * phi_b = a * phi_P + b 
 * c and d are the gradient coefficients that construct the boundary gradient
 * gradPhi_b = c * phi_P + d
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Project headers
#include "Integer.h"
#include "Scalar.h"

// ************************** struct BoundaryCoeffs ***************************

struct BoundaryCoeffs
{
    /// Construct coefficient arrays sized to numFaces (default: zeroGradient)
    explicit BoundaryCoeffs(Count numFaces = 0)
    :
        a(numFaces, S(1.0)),
        b(numFaces, S(0.0)),
        c(numFaces, S(0.0)),
        d(numFaces, S(0.0))
    {}

    ScalarList a; ///< Value internal
    ScalarList b; ///< Value boundary
    ScalarList c; ///< Gradient internal
    ScalarList d; ///< Gradient boundary
};
