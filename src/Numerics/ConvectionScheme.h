#ifndef CONVECTIONSCHEME_H
#define CONVECTIONSCHEME_H

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "CellData.h"

// Base class for all convection discretization schemes
class ConvectionDiscretization {
public:
    virtual ~ConvectionDiscretization() = default;
    virtual void getFluxCoefficients(
        Scalar F, // Convective mass flux rate through the face
        Scalar& a_P_conv,
        Scalar& a_N_conv,
        const Face& face,
        const std::vector<Cell>& cells
    ) const = 0;
};

// First-order Upwind Differencing Scheme (UDS)
class UpwindScheme : public ConvectionDiscretization {
public:
    void getFluxCoefficients(
        Scalar F,
        Scalar& a_P_conv,
        Scalar& a_N_conv,
        const Face& face,
        const std::vector<Cell>& cells
    ) const override;
};

// Central Differencing Scheme (CDS)
class CentralDifferenceScheme : public ConvectionDiscretization {
public:
    void getFluxCoefficients(
        Scalar F,
        Scalar& a_P_conv,
        Scalar& a_N_conv,
        const Face& face,
        const std::vector<Cell>& cells
    ) const override;
};

// Second-order Upwind Scheme
class SecondOrderUpwindScheme : public ConvectionDiscretization {
public:
    void getFluxCoefficients(
        Scalar F,
        Scalar& a_P_conv,
        Scalar& a_N_conv,
        const Face& face,
        const std::vector<Cell>& cells
    ) const override;

    // Calculate the second-order correction term
    Scalar calculateSecondOrderCorrection(
        const Face& face,
        const std::vector<Cell>& cells,
        const ScalarField& phi,
        Scalar F
    ) const;
};

#endif // CONVECTIONSCHEME_H