#ifndef CONVECTIONSCHEME_H
#define CONVECTIONSCHEME_H

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "CellData.h"
#include "FaceData.h"
#include "BoundaryConditions.h"
#include "BoundaryData.h"


// Base class for all convection discretization schemes
class ConvectionScheme {
public:
    virtual ~ConvectionScheme() = default;
    virtual void getFluxCoefficients(
        Scalar F, // Convective mass flux rate through the face
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const = 0;
};

/* First-order Upwind Differencing Scheme (UDS)
 * φ_f = φ_upwind
 * For F >= 0: φ_f = φ_P, so a_P = F, a_N = 0
 * For F < 0:  φ_f = φ_N, so a_P = 0, a_N = F
 */
class UpwindScheme : public ConvectionScheme {
public:
    void getFluxCoefficients(
        Scalar F,
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const override;
};

/* Bounded Central Differencing Scheme (BCDS) with gradient reformulation
 * φ_f = φ_P*w + φ_N*(1-w) + (∇φ_f · d_Pf)
 * Matrix uses upwind coefficients for stability
 * Higher-order accuracy via explicit correction term (deferred-correction)
 */
class CentralDifferenceScheme : public ConvectionScheme {
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
        Scalar F,
        const BoundaryConditions* bcManager = nullptr,
        const std::string& fieldName = ""
    ) const;

    // Calculate face value using gradient-based interpolation
    Scalar calculateFaceValue(
        const Face& face,
        const std::vector<Cell>& cells,
        const ScalarField& phi,
        const FaceVectorField& grad_phi_f,
        const BoundaryConditions* bcManager = nullptr,
        const std::string& fieldName = ""
    ) const;
};

/* Second-order Upwind Scheme with gradient-based formulation
 * φ_f = φ_P + (∇φ_P · d_Pf)
 * Matrix uses upwind coefficients for stability
 * Higher-order accuracy via explicit correction term (deferred-correction)
 */
class SecondOrderUpwindScheme : public ConvectionScheme {
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
        Scalar F,
        const BoundaryConditions* bcManager = nullptr,
        const std::string& fieldName = ""
    ) const;

    // Calculate face value using gradient-based upwind interpolation
    Scalar calculateFaceValue(
        const Face& face,
        const std::vector<Cell>& cells,
        const ScalarField& phi,
        const VectorField& grad_phi,
        Scalar F,
        const BoundaryConditions* bcManager = nullptr,
        const std::string& fieldName = ""
    ) const;
};

#endif