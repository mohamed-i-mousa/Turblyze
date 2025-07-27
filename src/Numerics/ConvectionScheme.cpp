#include "ConvectionScheme.h"
#include <algorithm>
#include "Face.h"
#include "Cell.h"

void UpwindScheme::getFluxCoefficients(
    Scalar F,
    Scalar& a_P_conv,
    Scalar& a_N_conv,
    const Face& face,
    const std::vector<Cell>& cells) const
{
    a_P_conv = std::max(F, 0.0);
    a_N_conv = std::min(F, 0.0);
}

void SecondOrderUpwindScheme::getFluxCoefficients(
    Scalar F,
    Scalar& a_P_conv,
    Scalar& a_N_conv,
    const Face& face,
    const std::vector<Cell>& cells) const
{
    // For the matrix assembly, use UDS for the implicit coefficients
    // This ensures diagonal dominance and stability
    a_P_conv = std::max(F, 0.0);
    a_N_conv = std::min(F, 0.0);
}

Scalar SecondOrderUpwindScheme::calculateSecondOrderCorrection(
    const Face& face,
    const std::vector<Cell>& cells,
    const ScalarField& phi,
    Scalar F) const
{
    if (face.isBoundary()) return S(0.0);

    size_t P = face.ownerCell;
    size_t N = face.neighbourCell.value();

    // Find the upwind cell and its upwind neighbor
    size_t U = (F >= S(0.0)) ? P : N;  // Upwind cell
    size_t D = (F >= S(0.0)) ? N : P;  // Downwind cell

    // Get the upwind cell's upwind neighbor
    // This would typically require additional mesh connectivity information
    // For now, we'll use a simplified approach
    const Cell& cell_U = cells[U];
    const Cell& cell_D = cells[D];

    // Calculate the gradient at the upwind cell
    Vector grad_phi_U = (phi[D] - phi[U]) * (cell_D.centroid - cell_U.centroid) / 
                       (cell_D.centroid - cell_U.centroid).magnitude();

    // Calculate the vector from upwind cell to face
    Vector d_Uf = face.centroid - cell_U.centroid;

    // Calculate the second-order correction
    // phi_f = phi_U + grad_phi_U . d_Uf
    Scalar phi_f_SOU = phi[U] + dot(grad_phi_U, d_Uf);

    // Calculate the first-order value (UDS)
    Scalar phi_f_UDS = phi[U];

    // Return the difference as the correction
    return F * (phi_f_SOU - phi_f_UDS);
}

void CentralDifferenceScheme::getFluxCoefficients(
    Scalar F,
    Scalar& a_P_conv,
    Scalar& a_N_conv,
    const Face& face,
    const std::vector<Cell>& cells) const
{
    if (face.isBoundary()) {
        a_P_conv = 0.0;  // Contribution to the owner cell
        a_N_conv = 0.0;  // Contribution to the virtual neighbour cell
        return;
    }

    // Get cell centroids
    const Vector& P = cells[face.ownerCell].centroid;
    const Vector& N = cells[face.neighbourCell.value()].centroid;
    const Vector& f = face.centroid;

    // Calculate weighting factor based on distance
    Scalar d_PN = (N - P).magnitude();
    Scalar d_Pf = (f - P).magnitude();

    Scalar w = d_Pf / d_PN; // Weight for cell N

    a_P_conv = F * (1.0 - w);
    a_N_conv = F * w;
}