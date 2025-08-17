#ifndef GRADIENTSCHEME_H
#define GRADIENTSCHEME_H

#include <vector>
#include <string>

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "FaceData.h"

/**
 * @file GradientScheme.h
 * @brief Gradient reconstruction and interpolation schemes
 * 
 * For each cell P, we solve the overdetermined system:
 *   phi_N - phi_P = ∇φ_P · (r_N - r_P) + O(|r_N - r_P|²)
 * 
 * This leads to the normal equations:  A^T A ∇φ_P = A^T b
 * where:
 *   A_ij = w_i * (r_i - r_P)_j
 *   b_i = w_i * (phi_i - phi_P)
 *   w_i = 1/|r_i - r_P|² (inverse distance squared weighting)
 *
 * @example
 * // 1. Calculate cell-centered gradients using least-squares
 * VectorField grad_phi = gradientScheme.LeastSquares(phi, allCells);
 * 
 * // 2. Interpolate gradients to face values
 * FaceVectorField grad_phi_faces = gradientScheme.interpolateGradientsToFaces(
 *     grad_phi, phi, allCells, allFaces, boundaryConditions, "phi");
 * 
 * // 3. Use the face gradients for flux calculations, etc.
 * for (size_t faceId = 0; faceId < allFaces.size(); ++faceId) {
 *     Vector grad_at_face = grad_phi_faces[faceId];
 *     // ... use grad_at_face for calculations
 * }
 */

/**
 * @brief Gradient reconstruction and interpolation schemes
 * 
 * This class provides methods for calculating cell-centered gradients
 * using least-squares reconstruction and interpolating them to face
 * values for use in flux calculations.
 */
class GradientScheme
{
public:
    GradientScheme() = default;

    /**
     * @brief Calculate cell-centered gradients using least-squares
     * @param phi Scalar field for gradient calculation
     * @param allCells Vector of all cells in the mesh
     * @return Vector field containing gradients at cell centers
     */
    VectorField LeastSquares
    (
        const ScalarField& phi,
        const std::vector<Cell>& allCells
    ) const;

    /**
     * @brief Interpolate cell-centered gradients to face values
     * * Interpolate cell-centered gradients to face values:
     * ∇φ_f = ∇φ_avg + [ (φ_N - φ_P)/|d_PN| - (∇φ_avg · e_PN) ] e_PN
     *
     * where:
     *   ∇φ_avg = g_P ∇φ_P + g_N ∇φ_N (distance-weighted average)
     *   e_PN   = d_PN / |d_PN| (unit vector P→N)
     *   g_P, g_N are interpolation weights (distance-based)
     *
     * This scheme preserves second-order accuracy.
     * 
     * @param grad_phi Cell-centered gradient field
     * @param phi Cell-centered scalar field
     * @param allCells Vector of all cells in the mesh
     * @param allFaces Vector of all faces in the mesh
     * @param boundaryConditions Boundary conditions manager
     * @param fieldName Name of the field for boundary condition lookup
     * @return Vector field containing gradients at face centers
     */
    FaceVectorField interpolateGradientsToFaces
    (
        const VectorField& grad_phi,
        const ScalarField& phi,
        const std::vector<Cell>& allCells,
        const std::vector<Face>& allFaces,
        const BoundaryConditions& boundaryConditions,
        const std::string& fieldName
    ) const;

private:
    /**
     * @brief Calculate boundary face gradient based on boundary condition type
     * @param face Boundary face
     * @param cellGradient Gradient at the owner cell
     * @param phi Scalar field values
     * @param boundaryConditions Boundary conditions manager
     * @param fieldName Name of the field for boundary condition lookup
     * @return Gradient vector at the boundary face
     */
    Vector calculateBoundaryFaceGradient
    (
        const Face& face,
        const Vector& cellGradient,
        const ScalarField& phi,
        const BoundaryConditions& boundaryConditions,
        const std::string& fieldName
    ) const;
};


#endif