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

// Rhie-Chow interpolation scheme for colocated grid momentum solution
class RhieChowScheme {
public:
    // Calculate face velocity using Rhie-Chow interpolation
    Vector calculateFaceVelocity(
        const Face& face,
        const std::vector<Cell>& cells,
        const ScalarField& u,
        const ScalarField& v,
        const ScalarField& w,
        const ScalarField& p,
        const VectorField& grad_u,
        const VectorField& grad_v,
        const VectorField& grad_w,
        const VectorField& grad_p,
        const ScalarField& a_u,
        const ScalarField& a_v,
        const ScalarField& a_w
    ) const;

    // Calculate pressure gradient correction term
    Vector calculatePressureGradientCorrection(
        const Face& face,
        const std::vector<Cell>& cells,
        const ScalarField& p,
        const VectorField& grad_p,
        const ScalarField& a_p
    ) const;
};

#endif // CONVECTIONSCHEME_H