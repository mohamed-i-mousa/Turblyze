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