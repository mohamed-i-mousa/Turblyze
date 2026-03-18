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
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const BoundaryConditions& bc
) noexcept
  : allFaces_(faces),
    allCells_(cells),
    bcManager_(bc)
{}


// ****************************** Public Methods ******************************

Vector GradientScheme::cellGradient
(
    const std::string& fieldName,
    const ScalarField& phi,
    size_t cellIndex,
    const FaceData<Scalar>* boundaryFaceValues
) const
{
    const Cell& cell = allCells_[cellIndex];

    Eigen::Matrix<Scalar,3,3> ATA;
    Eigen::Matrix<Scalar,3,1> ATb;
    Eigen::Matrix<Scalar,3,1> rVector;

    ATA.setZero();
    ATb.setZero();

    // Part 1: Internal neighbor cells contribution
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

        rVector << r.x(), r.y(), r.z();
        ATA.noalias() += w * (rVector * rVector.transpose());

        Scalar deltaPhi = phi[neighborIdx] - phi[cellIndex];
        ATb.noalias() += w * deltaPhi * rVector;
    }

    // Part 2: Boundary faces contribution
    for (size_t faceIdx : cell.faceIndices())
    {
        const Face& face = allFaces_[faceIdx];

        if (!face.isBoundary()) continue;  // Skip internal faces

        Vector r = face.centroid() - cell.centroid();
        Scalar rMagSqr = r.magnitudeSquared();

        Scalar w = S(1.0) / (rMagSqr + vSmallValue);

        rVector << r.x(), r.y(), r.z();
        ATA.noalias() += w * (rVector * rVector.transpose());

        Scalar phiBoundary = S(0.0);
        bool useOverride =
            boundaryFaceValues
         && face.idx() < boundaryFaceValues->size()
         && std::isfinite((*boundaryFaceValues)[face.idx()]);

        if (useOverride)
        {
            phiBoundary = (*boundaryFaceValues)[face.idx()];
        }
        else
        {
            // Get boundary value from BC
            phiBoundary =
                bcManager_.calculateBoundaryFaceValue
                (
                    fieldName,
                    phi,
                    face
                );
        }

        Scalar deltaPhi = phiBoundary - phi[cellIndex];
        ATb.noalias() += w * deltaPhi * rVector;
    }

    Eigen::LLT<Eigen::Matrix<Scalar,3,3>> llt(ATA);

    if (llt.info() == Eigen::Success)
    {
        Eigen::Matrix<Scalar,3,1> g = llt.solve(ATb);
        return Vector(g(0), g(1), g(2));
    }
    else
    {
        std::cerr
            << "WARNING: LLT failed for cell "
            << cellIndex
            << ", falling back to FullPivLU"
            << std::endl;

        Eigen::FullPivLU<Eigen::Matrix<Scalar,3,3>> lu(ATA);

        if (lu.isInvertible())
        {
            Eigen::Matrix<Scalar,3,1> g = lu.solve(ATb);
            return Vector(g(0), g(1), g(2));
        }
        else
        {
            throw
                std::runtime_error
                (
                    "Gradient computation failed"
                  + std::string(" for cell ")
                  + std::to_string(cellIndex)
                );
        }
    }
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
    size_t faceIndex
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
                face
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
    const Face& face
) const
{
    const BoundaryPatch* patch = face.patch();

    const BoundaryData* bc =
        bcManager_.fieldBC(patch->patchName(), fieldName);

    if (!bc) return cellGradient;

    switch (bc->type())
    {
        case BCType::NO_SLIP:
        case BCType::FIXED_VALUE:
        {
            Scalar boundaryValue = S(0.0);  // Default for NO_SLIP

            if (bc->type() == BCType::FIXED_VALUE)
            {
                if (bc->valueType() == BCValueType::SCALAR)
                {
                    boundaryValue = bc->fixedScalarValue();
                }
                else if (bc->valueType() == BCValueType::VECTOR)
                {
                    if (fieldName == "Ux")
                    {
                        boundaryValue = bc->vectorValue().x();
                    }
                    else if (fieldName == "Uy")
                    {
                        boundaryValue = bc->vectorValue().y();
                    }
                    else if (fieldName == "Uz")
                    {
                        boundaryValue = bc->vectorValue().z();
                    }
                    else
                    {
                        return cellGradient;
                    }
                }
                else
                {
                    return cellGradient;
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

        case BCType::K_WALL_FUNCTION:
        case BCType::NUT_WALL_FUNCTION:
        case BCType::OMEGA_WALL_FUNCTION:
        case BCType::ZERO_GRADIENT:
        {
            // Zero normal gradient: retain only tangential
            Vector tangentialGradient =
                cellGradient
              - dot(cellGradient, face.normal())
              * face.normal();

            return tangentialGradient;
        }

        case BCType::FIXED_GRADIENT:
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
