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

/*
 * Usage Example:
 * 
 * // 1. Calculate cell-centered gradients using least-squares
 * VectorField grad_phi = gradientScheme.LeastSquares(phi, allCells);
 * 
 * // 2. Interpolate gradients to face values
 * FaceVectorField grad_phi_faces = gradientScheme.interpolateGradientsToFaces(
 *     grad_phi, phi, allCells, allFaces);
 * 
 * // 3. Use the face gradients for flux calculations, non-orthogonal corrections, etc.
 * for (size_t faceId = 0; faceId < allFaces.size(); ++faceId) {
 *     Vector grad_at_face = grad_phi_faces[faceId];
 *     // ... use grad_at_face for calculations
 * }
 */

class GradientScheme {
public:
    GradientScheme() = default;

    // Cell-centred Least-Squares gradient reconstruction
    VectorField LeastSquares(
        const ScalarField& phi,
        const std::vector<Cell>& allCells
    ) const;

    // Interpolate cell-centered gradients to face values
    FaceVectorField interpolateGradientsToFaces(
        const VectorField& grad_phi,
        const ScalarField& phi,
        const std::vector<Cell>& allCells,
        const std::vector<Face>& allFaces
    ) const;
};


#endif