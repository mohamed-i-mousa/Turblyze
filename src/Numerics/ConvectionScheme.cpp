#include <algorithm>
#include <cmath>

#include "ConvectionScheme.h"
#include "GradientScheme.h"


// First-order Upwind Differencing Scheme (UDS)
void UpwindScheme::getFluxCoefficients(
    Scalar massFlowRate,
    Scalar& a_P_conv,
    Scalar& a_N_conv) const
{
    a_P_conv = std::max(massFlowRate, 0.0);
    a_N_conv = std::min(massFlowRate, 0.0);
}


// Bounded Central Differencing Scheme (BCDS) with gradient reformulation
void CentralDifferenceScheme::getFluxCoefficients(
    Scalar massFlowRate,
    Scalar& a_P_conv,
    Scalar& a_N_conv) const
{
    a_P_conv = std::max(massFlowRate, 0.0);
    a_N_conv = std::min(massFlowRate, 0.0);
}

Scalar CentralDifferenceScheme::calculateCentralDifferenceCorrection(
    const Face& face,
    const ScalarField& phi,
    const FaceVectorField& grad_phi_f,
    Scalar massFlowRate,
    const BoundaryConditions* bcManager,
    const std::string& fieldName) const
{
    if (face.isBoundary()) return 0.0;

    Scalar phi_face_central = calculateFaceValue(face, phi, grad_phi_f, bcManager, fieldName);

    size_t upwind_cell = (massFlowRate >= 0.0) ? face.ownerCell : face.neighbourCell.value();

    Scalar phi_face_UDS = phi[upwind_cell];
    
    // Return correction term: mdot * (φ_central - φ_upwind)
    return massFlowRate * (phi_face_central - phi_face_UDS);
}

Scalar CentralDifferenceScheme::calculateFaceValue(
    const Face& face,
    const ScalarField& phi,
    const FaceVectorField& grad_phi_f,
    const BoundaryConditions* bcManager,
    const std::string& fieldName) const
{
    if (face.isBoundary()) {
        if (bcManager && !fieldName.empty()) {
            return bcManager->calculateBoundaryFaceValue(face, phi, fieldName);
        }
        // Fallback to zero gradient if no BC info available
        std::cerr << "No BC specified for boundary face " << face.id << ". Defaulting to zero-gradient." << std::endl;
        return phi[face.ownerCell];
    }

    size_t P = face.ownerCell;
    size_t N = face.neighbourCell.value();
    
    // Calculate interpolation weight
    const Vector& d_Pf = face.d_Pf;
    Scalar d_P = face.d_Pf_mag;
    Scalar d_N = face.d_Nf_mag.value();
    Scalar w = d_N / (d_P + d_N);
    
    // Central difference formula: φ_f = φ_P*w + φ_N*(1-w) + (∇φ_f · d_Pf)
    Scalar phi_f = phi[P] * w + phi[N] * (S(1.0) - w) + dot(grad_phi_f[face.id], d_Pf);
    
    return phi_f;  
}

// Second-order Upwind Scheme with gradient-based formulation
void SecondOrderUpwindScheme::getFluxCoefficients(
    Scalar massFlowRate,
    Scalar& a_P_conv,
    Scalar& a_N_conv) const
{
    a_P_conv = std::max(massFlowRate, 0.0);
    a_N_conv = std::min(massFlowRate, 0.0);
}

Scalar SecondOrderUpwindScheme::calculateSecondOrderCorrection(
    const Face& face,
    const ScalarField& phi,
    const VectorField& grad_phi,
    Scalar massFlowRate,
    const BoundaryConditions* bcManager,
    const std::string& fieldName) const
{
    if (face.isBoundary()) return 0.0;

    Scalar phi_face_SOU = calculateFaceValue(face, phi, grad_phi, massFlowRate, bcManager, fieldName);
    
    size_t upwind_cell = (massFlowRate >= 0.0) ? face.ownerCell : face.neighbourCell.value();
    Scalar phi_face_UDS = phi[upwind_cell];
    
    // Return correction term: mdot * (φ_SOU - φ_UDS)
    return massFlowRate * (phi_face_SOU - phi_face_UDS);
}

Scalar SecondOrderUpwindScheme::calculateFaceValue(
    const Face& face,
    const ScalarField& phi,
    const VectorField& grad_phi,
    Scalar massFlowRate,
    const BoundaryConditions* bcManager,
    const std::string& fieldName) const
{
    if (face.isBoundary()) {
        if (bcManager && !fieldName.empty()) {
            return bcManager->calculateBoundaryFaceValue(face, phi, fieldName);
        } else {
            // Fallback to zero gradient if no BC info available
            std::cerr << "No BC specified for boundary face " << face.id << ". Defaulting to zero-gradient." << std::endl;
            return phi[face.ownerCell];
        }
    }

    // Determine upwind cell based on flow direction
    size_t upwind_cell = (massFlowRate >= 0.0) ? face.ownerCell : face.neighbourCell.value();
    
    // Second-order upwind formula: φ_f = φ_P + (∇φ_P · d_Pf)
    const Vector& d_Pf = face.d_Pf;
    Scalar phi_f = phi[upwind_cell] + dot(grad_phi[upwind_cell], d_Pf);
    
    return phi_f;
}