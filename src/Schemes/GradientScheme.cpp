/******************************************************************************
 * @file GradientScheme.cpp
 * @brief Implementation of gradient reconstruction schemes
 *****************************************************************************/

#include "GradientScheme.h"

#include <iostream>
#include <cmath>
#include <algorithm>

#include <omp.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>

#include "ErrorHandler.h"


// ************************* Special Member Functions *************************

GradientScheme::GradientScheme
(
    const Mesh& mesh,
    const BoundaryConditions& bc
)
:
    mesh_(mesh),
    bcManager_(bc)
{
    precomputeInverseATA();
}

// ****************************** Private Methods ******************************

void GradientScheme::precomputeInverseATA()
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

            const Scalar rMagSqr = r.magnitudeSquared();
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
            const Scalar rMagSqr = r.magnitudeSquared();
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

// ****************************** Public Methods ******************************

Vector GradientScheme::cellGradient
(
    const std::string& fieldName,
    const ScalarField& phi,
    size_t cellIndex,
    const FaceData<Scalar>* boundaryFaceValues
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

        const Scalar rMagSqr = r.magnitudeSquared();
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
        const Scalar rMagSqr = r.magnitudeSquared();
        const Scalar w = S(1.0) / (rMagSqr + smallValue);

        Scalar phiBoundary = S(0.0);
        const bool useOverride =
            boundaryFaceValues
         && face.idx() < boundaryFaceValues->size()
         && std::isfinite
            (
                (*boundaryFaceValues)[face.idx()]
            );

        if (useOverride)
        {
            phiBoundary =
                (*boundaryFaceValues)[face.idx()];
        }
        else
        {
            phiBoundary =
                bcManager_.boundaryFaceValue
                (
                    fieldName,
                    phi,
                    face
                );
        }

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


void GradientScheme::limitGradient
(
    const std::string& fieldName,
    const ScalarField& phi,
    VectorField& gradPhi
) const
{
    const size_t numCells = mesh_.numCells();
    
    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Cell& cell = mesh_.cells()[cellIdx];

        // Find min/max among cell P and all its neighbors
        Scalar phiMin = phi[cellIdx];
        Scalar phiMax = phi[cellIdx];

        for (size_t neighborIdx : cell.neighborCellIndices())
        {
            phiMin = std::min(phiMin, phi[neighborIdx]);
            phiMax = std::max(phiMax, phi[neighborIdx]);
        }

        // Include boundary face values so the limiter does not clip
        // physically correct near-boundary gradients (e.g. k near walls)
        for (size_t faceIdx : cell.faceIndices())
        {
            const Face& face = mesh_.faces()[faceIdx];
            if (!face.isBoundary()) continue;

            const Scalar phiBound =
                bcManager_.boundaryFaceValue
                (
                    fieldName,
                    phi,
                    face
                );
            phiMin = std::min(phiMin, phiBound);
            phiMax = std::max(phiMax, phiBound);
        }

        // Compute Barth-Jespersen limiter: min alpha over all faces
        Scalar alpha = S(1.0);

        for (size_t faceIdx : cell.faceIndices())
        {
            const Face& face = mesh_.faces()[faceIdx];
            const Vector r = face.centroid() - cell.centroid();
            const Scalar phiFace = phi[cellIdx] + dot(gradPhi[cellIdx], r);
            const Scalar delta = phiFace - phi[cellIdx];

            if (delta > smallValue)
            {
                alpha = std::min(alpha, (phiMax - phi[cellIdx]) / delta);
            }
            else if (delta < -smallValue)
            {
                alpha = std::min(alpha, (phiMin - phi[cellIdx]) / delta);
            }
        }

        alpha = std::clamp(alpha, S(0.0), S(1.0));

        gradPhi[cellIdx] = alpha * gradPhi[cellIdx];
    }
}


void GradientScheme::fieldGradient
(
    const std::string& fieldName,
    const ScalarField& phi,
    VectorField& gradPhi
) const
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradPhi[cellIdx] = cellGradient(fieldName, phi, cellIdx);
    }

    limitGradient(fieldName, phi, gradPhi);
}


Vector GradientScheme::faceGradient
(
    const std::string& fieldName,
    const ScalarField& phi,
    const Vector& gradPhiP,
    const Vector& gradPhiN,
    size_t faceIndex
) const
{
    const Face& face = mesh_.faces()[faceIndex];
    const size_t P = face.ownerCell();

    if (face.isBoundary())
    {
        return
            boundaryFaceGradient
            (
                fieldName,
                phi,
                gradPhiP,
                face
            );
    }
    else
    {
        const size_t N = face.neighborCell().value();
        const Vector dPN =
            mesh_.cells()[N].centroid() - mesh_.cells()[P].centroid();
        const Scalar dPNMag = dPN.magnitude();
        const Vector ePN = dPN / (dPNMag + vSmallValue);
        const Vector gradAvg = averageFaceGradient(face, gradPhiP, gradPhiN);
        const Scalar phiDiff = phi[N] - phi[P];
        const Scalar correction =
            (phiDiff / (dPNMag + vSmallValue))
          - dot(gradAvg, ePN);

        return gradAvg + correction * ePN;
    }
}


Vector GradientScheme::averageFaceGradient
(
    const Face& face,
    const Vector& gradPhiP,
    const Vector& gradPhiN
) const
{
    if (face.isBoundary())
    {
        FatalError
        (
            "GradientScheme::averageFaceGradient: "
            "Cannot average gradient on boundary face"
        );
    }

    const Scalar dPf = face.dPfMag();
    const Scalar dNf = face.dNfMag().value();
    const Scalar totalDist = dPf + dNf;

    const Scalar gP = dNf / (totalDist + vSmallValue);
    const Scalar gN = dPf / (totalDist + vSmallValue);

    return gP * gradPhiP + gN * gradPhiN;
}


Vector GradientScheme::boundaryFaceGradient
(
    const std::string& fieldName,
    const ScalarField& phi,
    const Vector& cellGradient,
    const Face& face
) const
{
    const BoundaryPatch& patch = face.patch()->get();

    const BoundaryData& bc =
        bcManager_.fieldBC(patch.patchName(), fieldName);

    using enum BCType;
    switch (bc.type())
    {
        case NO_SLIP:
        case FIXED_VALUE:
        {
            Scalar boundaryValue = S(0.0);  // Default for NO_SLIP

            if (bc.type() == FIXED_VALUE)
            {
                boundaryValue = bc.fixedScalarValue();
            }
            const Scalar cellValue = phi[face.ownerCell()];
            const Scalar dn = dot(face.dPf(), face.normal());
            const Scalar dPfMag = face.dPfMag();
            
            // Stabilization: clamp dn to minNormalFraction_ * ||dPf||
            const Scalar dnStabilized =
                std::max(dn, minNormalFraction_ * dPfMag);

            // ∂φ/∂n = (φ_boundary - φ_cell) / dnStabilized
            const Scalar normalGradient =
                (boundaryValue - cellValue)
              / dnStabilized;

            const Vector tangentialGradient =
                cellGradient
              - dot(cellGradient, face.normal())
              * face.normal();

            return tangentialGradient + normalGradient * face.normal();
        }

        case K_WALL_FUNCTION:
        case NUT_WALL_FUNCTION:
        case OMEGA_WALL_FUNCTION:
        case ZERO_GRADIENT:
        {
            // Zero normal gradient: retain only tangential
            const Vector tangentialGradient =
                cellGradient
              - dot(cellGradient, face.normal())
              * face.normal();

            return tangentialGradient;
        }

        case FIXED_GRADIENT:
        {
            const Scalar specifiedGradient = bc.fixedScalarGradient();

            // Project cell gradient onto tangential directions
            // and combine with specified normal gradient
            const Vector tangentialGradient =
                cellGradient
              - dot(cellGradient, face.normal())
              * face.normal();

            return tangentialGradient + specifiedGradient * face.normal();
        }

        default:
            return cellGradient;
    }
}
