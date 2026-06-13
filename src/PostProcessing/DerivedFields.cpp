/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file DerivedFields.cpp
 * @brief Implementation of derived cell-centered scalar fields
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "DerivedFields.h"

// Standard library headers
#include <cmath>

// Project headers
#include "Tensor.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

// ***************************** Internal Helpers *****************************

namespace
{

ScalarField computeMagnitude(const VectorField& field)
{
    ScalarField result;
    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < field.size(); ++cellIdx)
    {
        result[cellIdx] = magnitude(field[cellIdx]);
    }
    return result;
}

} // namespace


ScalarField velocityMagnitude
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz
)
{
    ScalarField result;
    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < Ux.size(); ++cellIdx)
    {
        result[cellIdx] = std::sqrt
        (
            Ux[cellIdx] * Ux[cellIdx]
          + Uy[cellIdx] * Uy[cellIdx]
          + Uz[cellIdx] * Uz[cellIdx]
        );
    }
    return result;
}

ScalarField vorticityMagnitude(const VectorField& vorticity)
{
    return computeMagnitude(vorticity);
}

ScalarField QCriterion
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
)
{
    ScalarField qCriterion;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < gradUx.size(); ++cellIdx)
    {
        // Q = 0.5 * (||Omega||^2 - ||S||^2)
        const Tensor gradU =
            tensorFromRows(gradUx[cellIdx], gradUy[cellIdx], gradUz[cellIdx]);

        const Scalar sMagSq = gradU.symm().magnitudeSquared();
        const Scalar oMagSq = gradU.skew().magnitudeSquared();

        qCriterion[cellIdx] = S(0.5) * (oMagSq - sMagSq);
    }

    return qCriterion;
}

ScalarField strainRateMagnitude
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
)
{
    ScalarField strainRateMag;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < gradUx.size(); ++cellIdx)
    {
        // Strain rate magnitude = sqrt(2 * S_ij * S_ij)
        const Tensor gradU = tensorFromRows
        (
            gradUx[cellIdx],
            gradUy[cellIdx],
            gradUz[cellIdx]
        );

        const Scalar symmMagSq = gradU.symm().magnitudeSquared();
        strainRateMag[cellIdx] = std::sqrt(S(2.0) * symmMagSq);
    }

    return strainRateMag;
}

} // namespace VTK
