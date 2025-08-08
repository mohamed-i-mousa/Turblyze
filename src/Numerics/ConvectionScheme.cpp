#include <algorithm>
#include <cmath>

#include "ConvectionScheme.h"
#include "GradientScheme.h"


// First-order Upwind Differencing Scheme (UDS)
void UpwindScheme::getFluxCoefficients(
    Scalar F,
    Scalar& a_P_conv,
    Scalar& a_N_conv) const
{
    a_P_conv = std::max(F, 0.0);
    a_N_conv = std::min(F, 0.0);
}


// Bounded Central Differencing Scheme (BCDS) with gradient reformulation
void CentralDifferenceScheme::getFluxCoefficients(
    Scalar F,
    Scalar& a_P_conv,
    Scalar& a_N_conv) const
{
    a_P_conv = std::max(F, 0.0);
    a_N_conv = std::min(F, 0.0);
}

Scalar CentralDifferenceScheme::calculateCentralDifferenceCorrection(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& phi,
    const FaceVectorField& grad_phi_f,
    Scalar F,
    const BoundaryConditions* bcManager,
    const std::string& fieldName) const
{
    if (face.isBoundary()) return 0.0;

    Scalar phi_face_central = calculateFaceValue(face, cells, phi, grad_phi_f, bcManager, fieldName);

    size_t upwind_cell = (F >= 0.0) ? face.ownerCell : face.neighbourCell.value();

    Scalar phi_face_upwind = phi[upwind_cell];
    
    // Return correction term: mdot * (φ_central - φ_upwind)
    return F * (phi_face_central - phi_face_upwind);
}

Scalar CentralDifferenceScheme::calculateFaceValue(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& phi,
    const FaceVectorField& grad_phi_f,
    const BoundaryConditions* bcManager,
    const std::string& fieldName) const
{
    if (face.isBoundary()) {
        if (bcManager && !fieldName.empty()) {
            return bcManager->calculateBoundaryFaceValue(face, phi, cells, fieldName);
        }
        // Fallback to zero gradient if no BC info available
        return phi[face.ownerCell];
    }

    size_t P = face.ownerCell;
    size_t N = face.neighbourCell.value();
    
    const Cell& cell_P = cells[P];
    const Cell& cell_N = cells[N];
    
    // Calculate interpolation weight
    Vector d_Nf = face.centroid - cell_N.centroid;
    Vector d_Pf = face.centroid - cell_P.centroid;
    Scalar d_P = d_Pf.magnitude();
    Scalar d_N = d_Nf.magnitude();
    Scalar w = d_N / (d_P + d_N);
    
    // Central difference formula: φ_f = φ_P*w + φ_N*(1-w) + (∇φ_f · d_Pf)
    Scalar phi_f = phi[P] * w + phi[N] * (S(1.0) - w) + dot(grad_phi_f[face.id], d_Pf);
    
    return phi_f;  
}

// Second-order Upwind Scheme with gradient-based formulation
void SecondOrderUpwindScheme::getFluxCoefficients(
    Scalar F,
    Scalar& a_P_conv,
    Scalar& a_N_conv) const
{
    a_P_conv = std::max(F, 0.0);
    a_N_conv = std::min(F, 0.0);
}

Scalar SecondOrderUpwindScheme::calculateSecondOrderCorrection(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& phi,
    const VectorField& grad_phi,
    Scalar F,
    const BoundaryConditions* bcManager,
    const std::string& fieldName) const
{
    if (face.isBoundary()) return 0.0;

    Scalar phi_face_SOU = calculateFaceValue(face, cells, phi, grad_phi, F, bcManager, fieldName);
    
    size_t upwind_cell = (F >= 0.0) ? face.ownerCell : face.neighbourCell.value();
    Scalar phi_face_UDS = phi[upwind_cell];
    
    // Return correction term: mdot * (φ_SOU - φ_UDS)
    return F * (phi_face_SOU - phi_face_UDS);
}

Scalar SecondOrderUpwindScheme::calculateFaceValue(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& phi,
    const VectorField& grad_phi,
    Scalar F,
    const BoundaryConditions* bcManager,
    const std::string& fieldName) const
{
    if (face.isBoundary()) {
        if (bcManager && !fieldName.empty()) {
            return bcManager->calculateBoundaryFaceValue(face, phi, cells, fieldName);
        } else {
            // Fallback to zero gradient if no BC info available
            return phi[face.ownerCell];
        }
    }

    // Determine upwind cell based on flow direction
    size_t upwind_cell = (F >= 0.0) ? face.ownerCell : face.neighbourCell.value();
    const Cell& cell_upwind = cells[upwind_cell];
    
    // Second-order upwind formula: φ_f = φ_P + (∇φ_P · d_Pf)
    Vector d_Pf = face.centroid - cell_upwind.centroid;
    Scalar phi_f = phi[upwind_cell] + dot(grad_phi[upwind_cell], d_Pf);
    
    return phi_f;
}