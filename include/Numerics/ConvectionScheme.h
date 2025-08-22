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


/**
 * @brief Base class for all convection discretization schemes
 * 
 * This class defines the interface for convection schemes used in
 * finite volume discretization of transport equations.
 */
class ConvectionScheme 
{
public:
    virtual ~ConvectionScheme() = default;
    
    /**
     * @brief Get convection flux coefficients
     * @param massFlowRate Convective mass flux rate through the face
     * @param a_P_conv Output coefficient for owner cell
     * @param a_N_conv Output coefficient for neighbor cell
     */
    virtual void getFluxCoefficients
    (
        Scalar massFlowRate,
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const = 0;
};

/**
 * @brief First-order Upwind Differencing Scheme (UDS)
 * 
 * φ_f = φ_upwind
 * For F >= 0: φ_f = φ_P, so a_P = F, a_N = 0
 * For F < 0: φ_f = φ_N, so a_P = 0, a_N = F
 */
class UpwindScheme : public ConvectionScheme {
public:
    /**
     * @brief Get upwind convection flux coefficients
     * @param massFlowRate Convective mass flux rate through the face
     * @param a_P_conv Output coefficient for owner cell
     * @param a_N_conv Output coefficient for neighbor cell
     */
    void getFluxCoefficients
    (
        Scalar massFlowRate,
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const override;
};

/**
 * @brief Bounded Central Differencing Scheme (BCDS) with gradient reformulation
 * 
 * φ_f = φ_P*w + φ_N*(1-w) + (∇φ_f · d_Pf)
 * Matrix uses upwind coefficients for stability
 * Higher-order accuracy via explicit correction term (deferred-correction)
 */
class CentralDifferenceScheme : public ConvectionScheme 
{
public:
    /**
     * @brief Get central difference convection flux coefficients
     * @param massFlowRate Convective mass flux rate through the face
     * @param a_P_conv Output coefficient for owner cell
     * @param a_N_conv Output coefficient for neighbor cell
     */
    void getFluxCoefficients
    (
        Scalar massFlowRate,
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const override;

    /**
     * @brief Calculate the correction term using gradients
     * @param face Face for interpolation
     * @param phi Cell-centered field values
     * @param grad_phi_f Face-centered gradients
     * @param massFlowRate Mass flow rate through face
     * @param bcManager Boundary conditions manager (optional)
     * @param fieldName Field name for boundary conditions (optional)
     * @return Correction term value
     */
    Scalar calculateCentralDifferenceCorrection
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& grad_phi_f,
        Scalar massFlowRate,
        const BoundaryConditions* bcManager = nullptr,
        const std::string& fieldName = ""
    ) const;

    /**
     * @brief Calculate face value using gradient-based interpolation
     * @param face Face for interpolation
     * @param phi Cell-centered field values
     * @param grad_phi_f Face-centered gradients
     * @param bcManager Boundary conditions manager (optional)
     * @param fieldName Field name for boundary conditions (optional)
     * @return Interpolated face value
     */
    Scalar calculateFaceValue
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& grad_phi_f,
        const BoundaryConditions* bcManager = nullptr,
        const std::string& fieldName = ""
    ) const;
};

/**
 * @brief Second-order Upwind Scheme with gradient-based formulation
 * 
 * φ_f = φ_P + (∇φ_P · d_Pf)
 * Matrix uses upwind coefficients for stability
 * Higher-order accuracy via explicit correction term (deferred-correction)
 */
class SecondOrderUpwindScheme : public ConvectionScheme
{
public:
    /**
     * @brief Get second-order upwind convection flux coefficients
     * @param massFlowRate Convective mass flux rate through the face
     * @param a_P_conv Output coefficient for owner cell
     * @param a_N_conv Output coefficient for neighbor cell
     */
    void getFluxCoefficients
    (
        Scalar massFlowRate,
        Scalar& a_P_conv,
        Scalar& a_N_conv
    ) const override;

    /**
     * @brief Calculate the second-order correction term using gradients
     * @param face Face for interpolation
     * @param phi Cell-centered field values
     * @param grad_phi Cell-centered gradient
     * @param massFlowRate Mass flow rate through face
     * @param bcManager Boundary conditions manager (optional)
     * @param fieldName Field name for boundary conditions (optional)
     * @return Correction term value
     */
    Scalar calculateSecondOrderCorrection
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& grad_phi_P,
        const Vector& grad_phi_N,        
        Scalar massFlowRate,
        const BoundaryConditions* bcManager = nullptr,
        const std::string& fieldName = ""
    ) const;

    /**
     * @brief Calculate face value using gradient-based upwind interpolation
     * @param face Face for interpolation
     * @param phi Cell-centered field values
     * @param grad_phi Cell-centered gradients
     * @param massFlowRate Mass flow rate through face
     * @param bcManager Boundary conditions manager (optional)
     * @param fieldName Field name for boundary conditions (optional)
     * @return Interpolated face value
     */
    Scalar calculateFaceValue
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& grad_phi_P,
        const Vector& grad_phi_N,        
        Scalar massFlowRate,
        const BoundaryConditions* bcManager = nullptr,
        const std::string& fieldName = ""
    ) const;
};

#endif