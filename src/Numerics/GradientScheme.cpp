/******************************************************************************
 * @file GradientScheme.cpp
 * @brief Implementation of gradient reconstruction schemes
 *****************************************************************************/

#include "GradientScheme.hpp"

#include <stdexcept>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>


// ****************************** Constructor ******************************

GradientScheme::GradientScheme
(
    std::span<const Face> faces,
    std::span<const Cell> cells,
    const BoundaryConditions& bc
)
  : allFaces_(faces),
    allCells_(cells),
    bcManager_(bc)
{
    precomputeInverseATA();
}


// ****************************** Private Methods ******************************

void GradientScheme::precomputeInverseATA()
{
    const size_t numCells = allCells_.size();
    invATA_.resize(numCells);

    Eigen::Matrix<Scalar,3,3> ATA;
    Eigen::Matrix<Scalar,3,1> rVector;

    size_t degenerateCells = 0;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Cell& cell = allCells_[cellIdx];

        ATA.setZero();

        // Neighbor cells contribution (purely geometric)
        for (size_t neighborIdx : cell.neighborCellIndices())
        {
            Vector r =
                allCells_[neighborIdx].centroid()
              - cell.centroid();

            Scalar rMagSqr = r.magnitudeSquared();
            Scalar w = S(1.0) / (rMagSqr + vSmallValue);

            rVector << r.x(), r.y(), r.z();
            ATA.noalias() += w * (rVector * rVector.transpose());
        }

        // Boundary faces contribution (purely geometric)
        for (size_t faceIdx : cell.faceIndices())
        {
            const Face& face = allFaces_[faceIdx];

            if (!face.isBoundary()) continue;

            Vector r = face.centroid() - cell.centroid();
            Scalar rMagSqr = r.magnitudeSquared();
            Scalar w = S(1.0) / (rMagSqr + vSmallValue);

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
        std::cerr
            << "WARNING: " << degenerateCells
            << " cells have degenerate least-squares"
            << " matrices (gradient will be zero)"
            << std::endl;
    }
}


// ****************************** Public Methods ******************************

Vector GradientScheme::cellGradient
(
    const std::string& fieldName,
    const ScalarField& phi,
    size_t cellIndex,
    const FaceData<Scalar>* boundaryFaceValues,
    std::optional<int> componentIdx
) const
{
    const Cell& cell = allCells_[cellIndex];

    Scalar b0 = S(0.0);
    Scalar b1 = S(0.0);
    Scalar b2 = S(0.0);

    // Part 1: Internal neighbor cells contribution to ATb
    for (size_t neighborIdx : cell.neighborCellIndices())
    {
        assert
        (
            neighborIdx < allCells_.size()
         && "Invalid neighbor ID - mesh topology corrupted"
        );

        const Cell& neighbor = allCells_[neighborIdx];
        Vector r = neighbor.centroid() - cell.centroid();

        Scalar rMagSqr = r.magnitudeSquared();
        Scalar w = S(1.0) / (rMagSqr + vSmallValue);

        Scalar wDeltaPhi =
            w * (phi[neighborIdx] - phi[cellIndex]);

        b0 += wDeltaPhi * r.x();
        b1 += wDeltaPhi * r.y();
        b2 += wDeltaPhi * r.z();
    }

    // Part 2: Boundary faces contribution to ATb
    for (size_t faceIdx : cell.faceIndices())
    {
        const Face& face = allFaces_[faceIdx];

        if (!face.isBoundary()) continue;

        Vector r = face.centroid() - cell.centroid();
        Scalar rMagSqr = r.magnitudeSquared();
        Scalar w = S(1.0) / (rMagSqr + vSmallValue);

        Scalar phiBoundary = S(0.0);
        bool useOverride =
            boundaryFaceValues
         && face.idx() < boundaryFaceValues->size()
         && std::isfinite(
                (*boundaryFaceValues)[face.idx()]);

        if (useOverride)
        {
            phiBoundary =
                (*boundaryFaceValues)[face.idx()];
        }
        else
        {
            phiBoundary =
                bcManager_.calculateBoundaryFaceValue
                (
                    fieldName,
                    phi,
                    face,
                    componentIdx
                );
        }

        Scalar wDeltaPhi =
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
    const ScalarField& phi,
    VectorField& gradPhi
) const noexcept
{
    for (size_t cellIdx = 0; cellIdx < allCells_.size(); ++cellIdx)
    {
        const Cell& cell = allCells_[cellIdx];
        Scalar phiP = phi[cellIdx];

        // Find min/max among cell P and all its neighbors
        Scalar phiMin = phiP;
        Scalar phiMax = phiP;

        for (size_t neighborIdx : cell.neighborCellIndices())
        {
            phiMin = std::min(phiMin, phi[neighborIdx]);
            phiMax = std::max(phiMax, phi[neighborIdx]);
        }

        // Compute Barth-Jespersen limiter: min alpha over all faces
        Scalar alpha = S(1.0);

        for (size_t faceIdx : cell.faceIndices())
        {
            const Face& face = allFaces_[faceIdx];
            Vector r = face.centroid() - cell.centroid();
            Scalar phiFace = phiP + dot(gradPhi[cellIdx], r);
            Scalar delta = phiFace - phiP;

            if (delta > vSmallValue)
            {
                alpha = std::min(alpha, (phiMax - phiP) / delta);
            }
            else if (delta < -vSmallValue)
            {
                alpha = std::min(alpha, (phiMin - phiP) / delta);
            }
        }

        alpha = std::clamp(alpha, S(0.0), S(1.0));

        gradPhi[cellIdx] = alpha * gradPhi[cellIdx];
    }
}

Vector GradientScheme::faceGradient
(
    const std::string& fieldName,
    const ScalarField& phi,
    const Vector& gradPhiP,
    const Vector& gradPhiN,
    size_t faceIndex,
    std::optional<int> componentIdx
) const
{
    const Face& face = allFaces_[faceIndex];
    size_t P = face.ownerCell();

    if (face.isBoundary())
    {
        return
            calculateBoundaryFaceGradient
            (
                fieldName,
                phi,
                gradPhiP,
                face,
                componentIdx
            );
    }
    else
    {
        size_t N = face.neighborCell().value();

        Vector dPN = allCells_[N].centroid() - allCells_[P].centroid();

        Scalar dPNMag = dPN.magnitude();

        Vector ePN = dPN / (dPNMag + vSmallValue);

        Vector gradAvg = averageFaceGradient(face, gradPhiP, gradPhiN);

        Scalar phiDiff = phi[N] - phi[P];

        Scalar correction =
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
        throw
            std::invalid_argument
            (
                "VectorLinearInterpolation: "
                "Cannot interpolate on boundary face"
            );
    }

    Scalar dPf = face.dPfMag();
    Scalar dNf = face.dNfMag().value();
    Scalar totalDist = dPf + dNf;

    Scalar gP = dNf / (totalDist + vSmallValue);
    Scalar gN = dPf / (totalDist + vSmallValue);

    return gP * gradPhiP + gN * gradPhiN;
}

Vector GradientScheme::calculateBoundaryFaceGradient
(
    const std::string& fieldName,
    const ScalarField& phi,
    const Vector& cellGradient,
    const Face& face,
    std::optional<int> componentIdx
) const
{
    const BoundaryPatch* patch = face.patch();

    const BoundaryData* bc =
        bcManager_.fieldBC(patch->patchName(), fieldName);

    if (!bc) return cellGradient;

    using enum BCType;
    switch (bc->type())
    {
        case NO_SLIP:
        case FIXED_VALUE:
        {
            Scalar boundaryValue = S(0.0);  // Default for NO_SLIP

            if (bc->type() == FIXED_VALUE)
            {
                if (componentIdx && bc->valueType() == BCValueType::VECTOR)
                {
                    const Vector& v = bc->vectorValue();
                    switch (*componentIdx)
                    {
                        case 0: boundaryValue = v.x(); break;
                        case 1: boundaryValue = v.y(); break;
                        case 2: boundaryValue = v.z(); break;
                    }
                }
                else if (bc->valueType() == BCValueType::SCALAR)
                {
                    boundaryValue = bc->fixedScalarValue();
                }
                else
                {
                    return cellGradient;
                }
            }
            else if (componentIdx && bc->valueType() == BCValueType::VECTOR)
            {
                // NO_SLIP with component index: extract from (0,0,0)
                const Vector& v = bc->vectorValue();
                switch (*componentIdx)
                {
                    case 0: boundaryValue = v.x(); break;
                    case 1: boundaryValue = v.y(); break;
                    case 2: boundaryValue = v.z(); break;
                }
            }

            Scalar cellValue = phi[face.ownerCell()];

            Scalar dn = dot(face.dPf(), face.normal());

            Scalar dPfMag = face.dPfMag();

            // Stabilization: limit dn to 5% minimum
            Scalar dnStabilized = std::max(dn, S(0.05) * dPfMag);

            // ∂φ/∂n = (φ_boundary - φ_cell) / dnStabilized
            Scalar normalGradient =
                (boundaryValue - cellValue)
              / dnStabilized;

            Vector tangentialGradient =
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
            Vector tangentialGradient =
                cellGradient
              - dot(cellGradient, face.normal())
              * face.normal();

            return tangentialGradient;
        }

        case FIXED_GRADIENT:
        {
            Scalar specifiedGradient = bc->fixedScalarGradient();

            // Project cell gradient onto tangential directions
            // and combine with specified normal gradient
            Vector tangentialGradient =
                cellGradient
              - dot(cellGradient, face.normal())
              * face.normal();

            return tangentialGradient + specifiedGradient * face.normal();
        }

        default:
            return cellGradient;
    }
}
