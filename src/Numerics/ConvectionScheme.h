#ifndef CONVECTIONSCHEME_H
#define CONVECTIONSCHEME_H

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "CellData.h"
#include "FaceData.h"

// Base class for all convection discretization schemes
class ConvectionDiscretization {
public:
    virtual ~ConvectionDiscretization() = default;
    virtual void getFluxCoefficients(
        Scalar F, // Convective mass flux rate through the face
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const = 0;
};

// First-order Upwind Differencing Scheme (UDS)
class UpwindScheme : public ConvectionDiscretization {
public:
    void getFluxCoefficients(
        Scalar F,
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const override;
};

// Bounded Central Differencing Scheme (BCDS) with gradient reformulation
class CentralDifferenceScheme : public ConvectionDiscretization {
public:
    void getFluxCoefficients(
        Scalar F,
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const override;

    // Calculate the correction term using gradients
    Scalar calculateCentralDifferenceCorrection(
        const Face& face,
        const std::vector<Cell>& cells,
        const ScalarField& phi,
        const FaceVectorField& grad_phi_f,

        Scalar F
    ) const;

    // Calculate face value using gradient-based interpolation
    Scalar calculateFaceValue(
        const Face& face,
        const std::vector<Cell>& cells,
        const ScalarField& phi,
        const FaceVectorField& grad_phi_f
    ) const;
};

// Second-order Upwind Scheme with gradient-based formulation
class SecondOrderUpwindScheme : public ConvectionDiscretization {
public:
    void getFluxCoefficients(
        Scalar F,
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const override;

    // Calculate the second-order correction term using gradients
    Scalar calculateSecondOrderCorrection(
        const Face& face,
        const std::vector<Cell>& cells,
        const ScalarField& phi,
        const VectorField& grad_phi,
        Scalar F
    ) const;

    // Calculate face value using gradient-based upwind interpolation
    Scalar calculateFaceValue(
        const Face& face,
        const std::vector<Cell>& cells,
        const ScalarField& phi,
        const VectorField& grad_phi,
        Scalar F
    ) const;
};

#endif // CONVECTIONSCHEME_H