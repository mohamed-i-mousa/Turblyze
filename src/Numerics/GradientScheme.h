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

    VectorField GreenGauss(
        const ScalarField& phi,
        const std::vector<Face>& allFaces,
        const std::vector<Cell>& allCells,
        const BoundaryConditions& bcManager
    ) const;

    // New: cell-centred Least-Squares gradient reconstruction
    VectorField LeastSquares(
        const ScalarField& phi,
        const std::vector<Cell>& allCells
    ) const;

private:
    // Helper to get the value of phi at a boundary face based on the BC
    Scalar getBoundaryFaceValue(
        const Face& face,
        size_t ownerCellId,
        const ScalarField& phi,
        const std::vector<Cell>& allCells,
        const BoundaryData& bc
    ) const;
};


#endif