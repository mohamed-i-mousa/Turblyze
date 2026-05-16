/******************************************************************************
 * @file DerivedFields.cpp
 * @brief Implementation of derived cell-centered scalar fields
 *****************************************************************************/

// Own header
#include "DerivedFields.h"

// STL includes
#include <cmath>

// Header includes
#include "Tensor.h"


namespace VTK
{

// ****************** Internal Helper Methods: Field Utilities *****************

ScalarField computeMagnitude(const VectorField& field)
{
    ScalarField magnitude;
    for (size_t idx = 0; idx < field.size(); ++idx)
    {
        magnitude[idx] = field[idx].magnitude();
    }
    return magnitude;
}

// *********************** Public API: Derived Fields **************************

ScalarField velocityMagnitude
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz
)
{
    ScalarField magnitude;
    for (size_t idx = 0; idx < Ux.size(); ++idx)
    {
        magnitude[idx] = std::sqrt
        (
            Ux[idx] * Ux[idx]
          + Uy[idx] * Uy[idx]
          + Uz[idx] * Uz[idx]
        );
    }
    return magnitude;
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

    for (size_t i = 0; i < gradUx.size(); ++i)
    {
        // Q = 0.5 * (||Omega||^2 - ||S||^2)
        const Tensor gradU = Tensor::fromRows(gradUx[i], gradUy[i], gradUz[i]);

        const Scalar sMagSq = gradU.symm().magnitudeSquared();
        const Scalar oMagSq = gradU.skew().magnitudeSquared();

        qCriterion[i] = S(0.5) * (oMagSq - sMagSq);
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

    for (size_t idx = 0; idx < gradUx.size(); ++idx)
    {
        // Strain rate magnitude = sqrt(2 * S_ij * S_ij)
        const Tensor gradU = Tensor::fromRows
        (
            gradUx[idx],
            gradUy[idx],
            gradUz[idx]
        );

        const Scalar symmMagSq = gradU.symm().magnitudeSquared();
        strainRateMag[idx] = std::sqrt(S(2.0) * symmMagSq);
    }

    return strainRateMag;
}

} // namespace VTK
