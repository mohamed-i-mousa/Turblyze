/******************************************************************************
 * @file LeastSquares.cpp
 * @brief Implementation of the weighted least-squares gradient scheme
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "LeastSquares.h"

// Standard library headers
#include <string>

// External library headers
#include <omp.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>

// Project headers
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
    size_t cellIndex
) const
{
    const Cell& cell = mesh_.cells()[cellIndex];

    Scalar b0 = S(0.0);
    Scalar b1 = S(0.0);
    Scalar b2 = S(0.0);

    // Part 1: Internal neighbor cells contribution to ATb
    for (size_t neighborIdx : cell.neighborCellIndices())
    {
        if (neighborIdx >= mesh_.numCells())
        {
            FatalError("Invalid neighbor ID - mesh topology corrupted");
        }

        const Cell& neighbor = mesh_.cells()[neighborIdx];
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
    for (size_t faceIdx : cell.faceIndices())
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (!face.isBoundary()) continue;

        const Vector r = face.centroid() - cell.centroid();
        const Scalar rMagSqr = magnitudeSquared(r);
        const Scalar w = S(1.0) / (rMagSqr + smallValue);

        const Scalar phiBoundary =
            bcManager_.boundaryFaceValue
            (
                field,
                phi,
                face
            );

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

// ****************************** Private Methods ******************************

void LeastSquares::precomputeInverseATA()
{
    const size_t numCells = mesh_.numCells();
    invATA_.resize(numCells);

    size_t degenerateCells = 0;

    #pragma omp parallel for schedule(static) reduction(+:degenerateCells)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Eigen::Matrix<Scalar,3,3> ATA;
        Eigen::Matrix<Scalar,3,1> rVector;
        const Cell& cell = mesh_.cells()[cellIdx];

        ATA.setZero();

        // Neighbor cells contribution (purely geometric)
        for (size_t neighborIdx : cell.neighborCellIndices())
        {
            const Vector r =
                mesh_.cells()[neighborIdx].centroid()
              - cell.centroid();

            const Scalar rMagSqr = magnitudeSquared(r);
            const Scalar w = S(1.0) / (rMagSqr + smallValue);

            rVector << r.x(), r.y(), r.z();
            ATA.noalias() += w * (rVector * rVector.transpose());
        }

        // Boundary faces contribution (purely geometric)
        for (size_t faceIdx : cell.faceIndices())
        {
            const Face& face = mesh_.faces()[faceIdx];

            if (!face.isBoundary()) continue;

            const Vector r = face.centroid() - cell.centroid();
            const Scalar rMagSqr = magnitudeSquared(r);
            const Scalar w = S(1.0) / (rMagSqr + smallValue);

            rVector << r.x(), r.y(), r.z();
            ATA.noalias() += w * (rVector * rVector.transpose());
        }

        // Invert ATA and store symmetric result
        Eigen::Matrix<Scalar,3,3> inv;
        bool inverted = false;

        Eigen::LLT<Eigen::Matrix<Scalar,3,3>> llt(ATA);

        if (llt.info() == Eigen::Success)
        {
            inv = llt.solve(Eigen::Matrix<Scalar,3,3>::Identity());
            inverted = true;
        }
        else
        {
            Eigen::FullPivLU<Eigen::Matrix<Scalar,3,3>>lu(ATA);

            if (lu.isInvertible())
            {
                inv = lu.inverse();
                inverted = true;
            }
        }

        if (inverted)
        {
            // Store upper triangle: {xx, xy, xz, yy, yz, zz}
            invATA_[cellIdx] =
            {
                inv(0,0), inv(0,1), inv(0,2),
                          inv(1,1), inv(1,2),
                                    inv(2,2)
            };
        }
        else
        {
            invATA_[cellIdx] = {0, 0, 0, 0, 0, 0};
            ++degenerateCells;
        }
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
