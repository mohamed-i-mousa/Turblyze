#include <algorithm>
#include <cmath>

#include "Face.h"
#include "Cell.h"
#include "ConvectionScheme.h"
#include "GradientScheme.h"

// ============================================================================
// FIRST-ORDER UPWIND SCHEME (UDS)
// ============================================================================
/*
 * First-order upwind: φ_f = φ_upwind
 * For F >= 0: φ_f = φ_P, so a_P = F, a_N = 0
 * For F < 0:  φ_f = φ_N, so a_P = 0, a_N = F
 */
void UpwindScheme::getFluxCoefficients(
    Scalar F,
    Scalar& a_P_conv,
    Scalar& a_N_conv) const
{
    a_P_conv = std::max(F, 0.0);
    a_N_conv = std::min(F, 0.0);
}

// ============================================================================
// CENTRAL DIFFERENCE SCHEME (CDS) - DEFERRED CORRECTION
// ============================================================================
/*
 * Central difference: φ_f = φ_P*w + φ_N*(1-w) + (∇φ_f · d_Pf)
 * Matrix uses upwind coefficients for stability
 * Higher-order accuracy via explicit correction term
 */
void CentralDifferenceScheme::getFluxCoefficients(
    Scalar F,
    Scalar& a_P_conv,
    Scalar& a_N_conv) const
{
    // Deferred correction: Use upwind coefficients for stability
    a_P_conv = std::max(F, 0.0);
    a_N_conv = std::min(F, 0.0);
}

Scalar CentralDifferenceScheme::calculateCentralDifferenceCorrection(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& phi,
    const FaceVectorField& grad_phi_f,
    Scalar F) const
{
    if (face.isBoundary()) return 0.0;

    // Calculate central difference face value
    Scalar phi_face_central = calculateFaceValue(face, cells, phi, grad_phi_f);
    // Calculate upwind face value
    size_t upwind_cell = (F >= 0.0) ? face.ownerCell : face.neighbourCell.value();
    Scalar phi_face_upwind = phi[upwind_cell];
    
    // Return correction term: F * (φ_central - φ_upwind)
    return F * (phi_face_central - phi_face_upwind);
}

Scalar CentralDifferenceScheme::calculateFaceValue(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& phi,
    const FaceVectorField& grad_phi_f) const
{
    if (face.isBoundary()) {
        return phi[face.ownerCell];
    }

    size_t P = face.ownerCell;
    size_t N = face.neighbourCell.value();
    
    const Cell& cell_P = cells[P];
    const Cell& cell_N = cells[N];
    
    // Calculate interpolation weight
    Vector d_PN = cell_N.centroid - cell_P.centroid;
    Vector d_Pf = face.centroid - cell_P.centroid;
    Scalar w = dot(d_Pf, d_PN) / dot(d_PN, d_PN);
    
    // Central difference formula: φ_f = φ_P*w + φ_N*(1-w) + (∇φ_f · d_Pf)
    Scalar phi_f = phi[P] * (S(1.0) - w) + phi[N] * w;
    
    // Use the provided face gradient instead of interpolating
    Scalar phi_correction = dot(grad_phi_f[face.id], d_Pf);
    
    return phi_f + phi_correction;  
}

// ============================================================================
// SECOND-ORDER UPWIND SCHEME (SOU) - DEFERRED CORRECTION
// ============================================================================
/*
 * Second-order upwind: φ_f = φ_P + (∇φ_P · d_Pf)
 * Matrix uses upwind coefficients for stability
 * Higher-order accuracy via explicit correction term
 */
void SecondOrderUpwindScheme::getFluxCoefficients(
    Scalar F,
    Scalar& a_P_conv,
    Scalar& a_N_conv) const
{
    // Use upwind coefficients for stability (deferred correction)
    a_P_conv = std::max(F, 0.0);
    a_N_conv = std::min(F, 0.0);
}

Scalar SecondOrderUpwindScheme::calculateSecondOrderCorrection(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& phi,
    const VectorField& grad_phi,
    Scalar F) const
{
    if (face.isBoundary()) return 0.0;

    // Calculate second-order upwind face value
    Scalar phi_face_SOU = calculateFaceValue(face, cells, phi, grad_phi, F);
    
    // Calculate first-order upwind face value
    size_t upwind_cell = (F >= 0.0) ? face.ownerCell : face.neighbourCell.value();
    Scalar phi_face_UDS = phi[upwind_cell];
    
    // Return correction term: F * (φ_SOU - φ_UDS)
    return F * (phi_face_SOU - phi_face_UDS);
}

Scalar SecondOrderUpwindScheme::calculateFaceValue(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& phi,
    const VectorField& grad_phi,
    Scalar F) const
{
    if (face.isBoundary()) {
        return phi[face.ownerCell];
    }

    // Determine upwind cell based on flow direction
    size_t upwind_cell = (F >= 0.0) ? face.ownerCell : face.neighbourCell.value();
    const Cell& cell_upwind = cells[upwind_cell];
    
    // Second-order upwind formula: φ_f = φ_P + (∇φ_P · d_Pf)
    Vector d_Pf = face.centroid - cell_upwind.centroid;
    Scalar phi_f = phi[upwind_cell] + dot(grad_phi[upwind_cell], d_Pf);
    
    return phi_f;
}

// ============================================================================
// RHIE-CHOW INTERPOLATION SCHEME
// ============================================================================

Vector RhieChowScheme::calculateFaceVelocity(
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
    const ScalarField& a_w) const
{
    if (face.isBoundary()) {
        // For boundary faces, use the owner cell velocity
        return Vector(u[face.ownerCell], v[face.ownerCell], w[face.ownerCell]);
    }

    size_t P = face.ownerCell;
    size_t N = face.neighbourCell.value();
    
    const Cell& cell_P = cells[P];
    const Cell& cell_N = cells[N];
    
    // Calculate interpolation weight
    Vector d_PN = cell_N.centroid - cell_P.centroid;
    Vector d_Pf = face.centroid - cell_P.centroid;
    Scalar weight = dot(d_Pf, d_PN) / dot(d_PN, d_PN);
    
    // Linear interpolation of velocity components
    Scalar u_face = (S(1.0) - weight) * u[P] + weight * u[N];
    Scalar v_face = (S(1.0) - weight) * v[P] + weight * v[N];
    Scalar w_face = (S(1.0) - weight) * w[P] + weight * w[N];
    
    // Linear interpolation of pressure gradients
    Vector grad_p_face = (S(1.0) - weight) * grad_p[P] + weight * grad_p[N];
    
    // Linear interpolation of coefficient fields
    Scalar a_u_face = (S(1.0) - weight) * a_u[P] + weight * a_u[N];
    Scalar a_v_face = (S(1.0) - weight) * a_v[P] + weight * a_v[N];
    Scalar a_w_face = (S(1.0) - weight) * a_w[P] + weight * a_w[N];
    
    // Calculate pressure gradient correction
    Vector pressure_correction = calculatePressureGradientCorrection(
        face, cells, p, grad_p, a_u); // Using a_u field
    
    // Apply Rhie-Chow correction
    Vector velocity_correction = pressure_correction / a_u_face;
    
    return Vector(u_face + velocity_correction.x, 
                  v_face + velocity_correction.y, 
                  w_face + velocity_correction.z);
}

Vector RhieChowScheme::calculatePressureGradientCorrection(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& p,
    const VectorField& grad_p,
    const ScalarField& a_p) const
{
    if (face.isBoundary()) return Vector(0.0, 0.0, 0.0);
    
    size_t P = face.ownerCell;
    size_t N = face.neighbourCell.value();
    
    const Cell& cell_P = cells[P];
    const Cell& cell_N = cells[N];
    
    // Calculate interpolation weight
    Vector d_PN = cell_N.centroid - cell_P.centroid;
    Vector d_Pf = face.centroid - cell_P.centroid;
    Scalar weight = dot(d_Pf, d_PN) / dot(d_PN, d_PN);
    
    // Interpolate pressure gradient to face
    Vector grad_p_face = (S(1.0) - weight) * grad_p[P] + weight * grad_p[N];
    
    // Calculate pressure difference across the face
    Vector face_normal = face.area * face.normal;
    Scalar p_face = (S(1.0) - weight) * p[P] + weight * p[N];
    
    // Rhie-Chow correction term
    // This is a simplified version - in practice, you'd need the full
    // momentum equation coefficients and more sophisticated interpolation
    Vector correction = Vector(0.0, 0.0, 0.0);
    
    // The correction should account for the difference between
    // interpolated pressure gradient and actual pressure difference
    // This is a placeholder for the full implementation
    Scalar a_p_face = (S(1.0) - weight) * a_p[P] + weight * a_p[N];
    
    return correction;
}