/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LeastSquares.cpp
 * @brief Implementation of the weighted least-squares gradient scheme
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "LeastSquares.h"

// Project headers
#include "Cramer3x3.h"
#include "ErrorHandler.h"

// ************************* Special Member Functions *************************

LeastSquares::LeastSquares
(
    const Mesh& mesh,
    const BoundaryConditions& bc
)
:
    GradientScheme(mesh, bc)
{
    precomputeInverseATA();
}

// ****************************** Public Methods ******************************

Vector LeastSquares::cellGradient
(
    Field field,
    const ScalarField& phi,
    Index cellIndex
) const
{
    const Cell& cell = mesh().cells()[cellIndex];

    Scalar b0 = S(0.0);
    Scalar b1 = S(0.0);
    Scalar b2 = S(0.0);

    // Part 1: Internal neighbor cells contribution to ATb
    for (Index neighborIdx : cell.neighborCellIndices())
    {
        const Cell& neighbor = mesh().cells()[neighborIdx];
        const Vector r = neighbor.centroid() - cell.centroid();

        const Scalar rMagSqr = magnitudeSquared(r);
        const Scalar w = S(1.0) / (rMagSqr + smallValue);

        const Scalar wDeltaPhi =
            w * (phi[neighborIdx] - phi[cellIndex]);

        b0 += wDeltaPhi * r.x();
        b1 += wDeltaPhi * r.y();
        b2 += wDeltaPhi * r.z();
    }

    // Part 2: Boundary faces contribution to ATb
    for (Index faceIdx : cell.faceIndices())
    {
        const Face& f = mesh().faces()[faceIdx];

        if (!f.isBoundary()) continue;

        const Vector r = f.centroid() - cell.centroid();
        const Scalar rMagSqr = magnitudeSquared(r);
        const Scalar w = S(1.0) / (rMagSqr + smallValue);

        const Scalar phiBoundary =
            bcManager().faceValue(f, phi[f.ownerCell()], field);

        const Scalar wDeltaPhi =
            w * (phiBoundary - phi[cellIndex]);

        b0 += wDeltaPhi * r.x();
        b1 += wDeltaPhi * r.y();
        b2 += wDeltaPhi * r.z();
    }

    // g = inv(ATA) * ATb  (symmetric 3x3 mat-vec multiply)
    const auto& inv = invATA_[cellIndex];

    return Vector
    (
        inv[0]*b0 + inv[1]*b1 + inv[2]*b2,
        inv[1]*b0 + inv[3]*b1 + inv[4]*b2,
        inv[2]*b0 + inv[4]*b1 + inv[5]*b2
    );
}


// ****************************** Private Methods *****************************

void LeastSquares::precomputeInverseATA()
{
    // Sized over every cell; ghosts are not gradient sites and stay zero
    invATA_.resize(mesh().numCells());

    const Count numDomainCells = mesh().numDomainCells();

    Count degenerateCells = 0;

    for (Index cellIdx = 0; cellIdx < numDomainCells; ++cellIdx)
    {
        const Cell& cell = mesh().cells()[cellIdx];

        // Assemble ATA in row-major flat array (symmetric 3x3)
        Scalar ATA[9] = {};

        // Neighbor cells contribution (purely geometric)
        for (Index neighborIdx : cell.neighborCellIndices())
        {
            const Vector r =
                mesh().cells()[neighborIdx].centroid()
              - cell.centroid();

            const Scalar rMagSqr = magnitudeSquared(r);
            const Scalar w = S(1.0) / (rMagSqr + smallValue);

            ATA[0] += w * r.x() * r.x();
            ATA[1] += w * r.x() * r.y();
            ATA[2] += w * r.x() * r.z();
            ATA[4] += w * r.y() * r.y();
            ATA[5] += w * r.y() * r.z();
            ATA[8] += w * r.z() * r.z();
        }

        // Boundary faces contribution (purely geometric)
        for (Index faceIdx : cell.faceIndices())
        {
            const Face& f = mesh().faces()[faceIdx];

            if (!f.isBoundary()) continue;

            const Vector r = f.centroid() - cell.centroid();
            const Scalar rMagSqr = magnitudeSquared(r);
            const Scalar w = S(1.0) / (rMagSqr + smallValue);

            ATA[0] += w * r.x() * r.x();
            ATA[1] += w * r.x() * r.y();
            ATA[2] += w * r.x() * r.z();
            ATA[4] += w * r.y() * r.y();
            ATA[5] += w * r.y() * r.z();
            ATA[8] += w * r.z() * r.z();
        }

        // Fill symmetric lower triangle
        ATA[3] = ATA[1];
        ATA[6] = ATA[2];
        ATA[7] = ATA[5];

        // Solve for each column of the inverse
        const Scalar e0[3] = {S(1.0), S(0.0), S(0.0)};
        const Scalar e1[3] = {S(0.0), S(1.0), S(0.0)};
        const Scalar e2[3] = {S(0.0), S(0.0), S(1.0)};

        Scalar col0[3], col1[3], col2[3];
        solve3x3(ATA, e0, col0);
        solve3x3(ATA, e1, col1);
        solve3x3(ATA, e2, col2);

        // Store upper triangle: {xx, xy, xz, yy, yz, zz}
        invATA_[cellIdx] =
        {
            col0[0], col1[0], col2[0],
                     col1[1], col2[1],
                              col2[2]
        };
    }

    if (degenerateCells > 0)
    {
        Warning
        (
            std::to_string(degenerateCells)
          + " cells have degenerate least-squares"
            " matrices (gradient will be zero)"
        );
    }
}
