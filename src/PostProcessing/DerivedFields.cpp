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

ScalarField velocityMagnitude(const VectorField& velocity)
{
    return computeMagnitude(velocity);
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
        Tensor gradU = Tensor::fromRows(gradUx[i], gradUy[i], gradUz[i]);

        Scalar sMagSq = gradU.symm().magnitudeSquared();
        Scalar oMagSq = gradU.skew().magnitudeSquared();

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
        Tensor gradU = Tensor::fromRows
        (
            gradUx[idx],
            gradUy[idx],
            gradUz[idx]
        );

        Scalar symmMagSq = gradU.symm().magnitudeSquared();
        strainRateMag[idx] = std::sqrt(S(2.0) * symmMagSq);
    }

    return strainRateMag;
}

} // namespace VTK
