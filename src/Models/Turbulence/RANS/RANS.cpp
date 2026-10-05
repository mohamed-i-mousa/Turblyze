/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file RANS.cpp
 * @brief Shared two-equation RANS services and residual bookkeeping
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "RANS.h"

// Standard library headers
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

// Project headers
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "Field.h"
#include "HaloExchange.h"
#include "Matrix.h"
#include "Reduce.h"
#include "LinearInterpolation.h"
#include "TimeScheme.h"

// ************************* Special Member Functions *************************

RANS::RANS
(
    const Mesh& mesh,
    const BoundaryConditions& bc,
    const TimeScheme& timeScheme,
    const GradientScheme& gradScheme,
    const ConvectionScheme& kScheme,
    LinearSolver& kSolver,
    const ConvectionScheme& dissipationScheme,
    LinearSolver& dissipationSolver,
    Scalar deltaT,
    Scalar nu,
    Scalar alphaK,
    Scalar alphaDissipation,
    bool debug
)
:
    mesh_{mesh},
    bcManager_{bc},
    timeScheme_{timeScheme},
    gradientScheme_{gradScheme},
    matrixConstruct_{std::make_unique<Matrix>(mesh_, bcManager_)},
    kConvectionScheme_{kScheme},
    kSolver_{kSolver},
    dissipationConvectionScheme_{dissipationScheme},
    dissipationSolver_{dissipationSolver},
    nu_{nu},
    deltaT_{deltaT},
    alphaK_{alphaK},
    alphaDissipation_{alphaDissipation},
    debug_{debug}
{}

RANS::~RANS() noexcept = default;

// ***************************** Transient Stepping **************************

void RANS::beginTimeStep()
{
    // Snapshot the converged previous-step fields as phi^n
    kPrevStep_ = k_;
    dissipationPrevStep_ = dissipation();
}


void RANS::updatePrevStepDerivatives()
{
    const ScalarField& dissipationNew = dissipation();
    const Count numCells = mesh_.numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar volume = mesh_.cells()[cellIdx].volume();

        kDdtPrevStep_[cellIdx] = timeScheme_.updateDdtPrevStep
        (
            volume,
            deltaT_,
            k_[cellIdx],
            kPrevStep_[cellIdx],
            kDdtPrevStep_[cellIdx]
        );
        dissipationDdtPrevStep_[cellIdx] = timeScheme_.updateDdtPrevStep
        (
            volume,
            deltaT_,
            dissipationNew[cellIdx],
            dissipationPrevStep_[cellIdx],
            dissipationDdtPrevStep_[cellIdx]
        );
    }
}


// ************************ Inlet Condition Calculators ***********************

Scalar RANS::inletK
(
    const Vector& velocity,
    Scalar turbulenceIntensity
) noexcept
{
    const Scalar uPrime = turbulenceIntensity * magnitude(velocity);
    return std::max(S(1.5) * uPrime * uPrime, smallValue);
}

// ***************************** Accessor Methods *****************************

Scalar RANS::boundaryTurbulentViscosity(const Face& face) const
{
    const BoundaryPatch& patch = *face.patch();

    if (bcManager_.boundaryType(patch.name(), Field::nut).isWallModelled())
    {
        return nutWall_[face.idx()];
    }

    return TurbulenceModel::boundaryTurbulentViscosity(face);
}


FaceData<Scalar> RANS::wallShearStress
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz
) const
{
    FaceData<Scalar> shearStress(mesh_, S(0.0));

    for (Index i = 0; i < wallFunctionFaceIndices_.size(); ++i)
    {
        const Index faceIdx = wallFunctionFaceIndices_[i];
        const auto& face = mesh_.faces()[faceIdx];
        const Index cellIdx = face.ownerCell();

        // Project velocity onto wall plane (tangential component)
        const Vector Ucell(Ux[cellIdx], Uy[cellIdx], Uz[cellIdx]);
        const Scalar normalVelocity = dot(Ucell, face.normal());
        const Vector tangentVelocity =
            Ucell - face.normal() * normalVelocity;
        const Scalar tangentVelocityMag = magnitude(tangentVelocity);

        const Scalar tau =
            yPlus_[face.idx()] < yPlusLam_
          ? nu_ * tangentVelocityMag / y_[face.idx()]
          : cmu25() * cmu25() * k_[cellIdx];

        shearStress[face.idx()] = std::min(tau, S(1000.0));
    }

    return shearStress;
}


RANS::CellDataPair RANS::cellDataOutputs() const
{
    return
    {
        {"k", &k_},
        {dissipationName(), &dissipation()},
        {"nut", &nut_},
        {"wallDistance", &wallDistance_}
    };
}


RANS::BoundaryDataPair RANS::boundaryDataOutputs() const
{
    return
    {
        {"yPlus", &yPlus_}
    };
}


RANS::ResidualPair RANS::residualOutputs() const
{
    return
    {
        {"k", lastKResidual_},
        {dissipationName(), lastDissipationResidual_}
    };
}

// ****************************** Shared Methods ******************************

std::optional<TransientTerm> RANS::ddtTerm
(
    const ScalarField& phiPrevStep,
    const ScalarField& ddtPrevStep
) const
{
    if (!timeScheme_.isTransient())
    {
        return std::nullopt;
    }

    return TransientTerm{timeScheme_, deltaT_, phiPrevStep, &ddtPrevStep};
}


void RANS::cellToFaceDiffusion
(
    const ScalarField& cellGamma,
    FaceFluxField& faceGamma
) const
{
    const Count numFaces = mesh_.numFaces();

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            faceGamma[faceIdx] = cellGamma[face.ownerCell()];
        }
        else
        {
            faceGamma[faceIdx] = interpolateToFace(mesh_, face, cellGamma);
        }
    }
}


void RANS::updateWallDistance()
{
    wallDistanceConverged_ = false;
    wallDistance_.setAll(S(1e10));
    nearestWallPoint_.setAll(Vector{S(1e15), S(1e15), S(1e15)});

    // Seed wall-adjacent cells with the distance to each wall centroid
    for (const auto& face : mesh_.faces())
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch& patch = *face.patch();

        if (patch.type() == PatchType::processor || !bcManager_.boundaryType(patch, Field::Ux).isWall())
        {
            continue;
        }

        const Index cellIdx = face.ownerCell();
        const Vector cellCenter = mesh_.cells()[cellIdx].centroid();
        const Vector faceCenter = face.centroid();
        const Vector normal = face.normal();
        const Vector cellToFace = faceCenter - cellCenter;
        const Scalar dist = std::abs(dot(cellToFace, normal));

        if (dist < wallDistance_[cellIdx])
        {
            wallDistance_[cellIdx] = dist;
            nearestWallPoint_[cellIdx] = faceCenter;
        }
    }

    // Iterative propagation
    const Count maxIterations = 100;
    const Scalar tolerance = smallValue;

    for (Index iter = 0; iter < maxIterations; ++iter)
    {
        // The wave crosses the cuts one exchange per sweep
        Halo::exchange({&wallDistance_});
        Halo::exchange({&nearestWallPoint_});

        Scalar maxChange = S(0.0);

        for (const auto& face : mesh_.faces())
        {
            if (face.isBoundary()) continue;

            const Index owner = face.ownerCell();
            const Index neighbor = face.neighborCell().value();
            const Vector ownerCenter = mesh_.cells()[owner].centroid();
            const Vector neighborCenter = mesh_.cells()[neighbor].centroid();

            // Try to improve owner using neighbor's nearest wall point
            {
                const Vector candidatePoint = nearestWallPoint_[neighbor];
                const Scalar newDist =
                    magnitude(ownerCenter - candidatePoint);

                if (newDist < wallDistance_[owner])
                {
                    const Scalar change = wallDistance_[owner] - newDist;
                    maxChange = std::max(maxChange, change);
                    wallDistance_[owner] = newDist;
                    nearestWallPoint_[owner] = candidatePoint;
                }
            }

            // Try to improve neighbor using owner's nearest wall point
            {
                const Vector candidatePoint = nearestWallPoint_[owner];
                const Scalar newDist =
                    magnitude(neighborCenter - candidatePoint);

                if (newDist < wallDistance_[neighbor])
                {
                    const Scalar change = wallDistance_[neighbor] - newDist;
                    maxChange = std::max(maxChange, change);
                    wallDistance_[neighbor] = newDist;
                    nearestWallPoint_[neighbor] = candidatePoint;
                }
            }
        }

        // Every rank sweeps the same number of times or this deadlocks
        if (globalMax(maxChange) < tolerance)
        {
            wallDistanceConverged_ = true;
            break;
        }
    }
}


void RANS::updateYPlusLam(Scalar kappa, Scalar E)
{
    Scalar yPlusLam = S(11.0);

    // 10 iterations to solve yPlusLam = ln(E*yPlusLam) / kappa
    for (int i = 0; i < 10; ++i)
    {
        yPlusLam =
            std::log(std::max(E * yPlusLam, S(1.0))) / kappa;
    }

    yPlusLam_ = yPlusLam;
}


void RANS::initializeWallFunctionGeometry
(
    const BoundaryConditions& bcManager,
    Field wallFunctionField
)
{
    const Count numCells = mesh_.numDomainCells();

    wallFunctionFaceIndices_.clear();
    wallCellIndices_.clear();
    wallCellFraction_.clear();
    wallFaceWeight_.setAll(S(0.0));
    y_.setAll(S(0.0));

    CountList wallFaceCountPerCell(numCells, 0);
    ScalarList totalWallArea(numCells, S(0.0));

    for (Index faceIdx = 0; faceIdx < mesh_.numFaces(); ++faceIdx)
    {
        const auto& face = mesh_.faces()[faceIdx];
        if (!face.isBoundary()) continue;

        const BoundaryPatch& patch = *face.patch();
        if (patch.type() == PatchType::processor || !bcManager.boundaryType(patch, Field::Ux).isWall())
        {
            continue;
        }

        const BoundaryType& bc =
            bcManager.boundaryType(patch.name(), wallFunctionField);
        if (!bc.isWallModelled()) continue;

        wallFunctionFaceIndices_.push_back(faceIdx);

        const Index cellIdx = face.ownerCell();
        totalWallArea[cellIdx] += face.projectedArea();
        ++wallFaceCountPerCell[cellIdx];
    }

    // Compute per-face weight = faceArea / totalWallArea
    for (Index faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const Index cellIdx = face.ownerCell();
        const Scalar area = totalWallArea[cellIdx];

        if (area > S(0.0))
        {
            wallFaceWeight_[face.idx()] = face.projectedArea() / area;
        }

        // Cache the owner-cell wall-normal distance
        y_[face.idx()] = std::max
        (
            std::abs(dot(mesh_.dPf(face), face.normal())),
            vSmallValue
        );
    }

    // Build unique wall cell indices
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        if (wallFaceCountPerCell[cellIdx] > 0)
        {
            wallCellIndices_.push_back(cellIdx);
        }
    }

    // Compute wallCellFraction = wallFunctionArea / totalPolyWallArea
    ScalarList totalPolyWallArea(numCells, S(0.0));

    for (Index cellIdx : wallCellIndices_)
    {
        for (Index faceIdx : mesh_.cells()[cellIdx].faceIndices())
        {
            const auto& face = mesh_.faces()[faceIdx];
            if (!face.isBoundary()) continue;
            
            if
            (
                face.patch()->type() != PatchType::processor
             && bcManager.boundaryType(*face.patch(), Field::Ux).isWall()
            )
            {
                totalPolyWallArea[cellIdx] += face.projectedArea();
            }
        }
    }

    wallCellFraction_.resize(wallCellIndices_.size());

    constexpr Scalar wallCellFractionTol = S(0.1);

    for (Index i = 0; i < wallCellIndices_.size(); ++i)
    {
        const Index cellIdx = wallCellIndices_[i];
        Scalar rawFraction = S(1.0);

        if (totalPolyWallArea[cellIdx] > S(0.0))
        {
            rawFraction =
                totalWallArea[cellIdx] / totalPolyWallArea[cellIdx];
        }

        wallCellFraction_[i] =
            std::max
            (
                (rawFraction - wallCellFractionTol)
              / (S(1.0) - wallCellFractionTol),
                S(0.0)
            );
    }
}


void RANS::updateYPlus()
{
    for (Index i = 0; i < wallFunctionFaceIndices_.size(); ++i)
    {
        const Index faceIdx = wallFunctionFaceIndices_[i];
        const auto& face = mesh_.faces()[faceIdx];
        const Index cellIdx = face.ownerCell();

        yPlus_[face.idx()] =
            cmu25() * std::sqrt(k_[cellIdx]) * y_[face.idx()] / nu_;
    }
}


ScalarField RANS::strainRateSquared(const TensorField& gradU) const
{
    const Count numCells = mesh_.numDomainCells();

    ScalarField strainRateSq{mesh_};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // S² = 2 S_ij S_ij, where S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        const Scalar symmMagSq = gradU[cellIdx].symm().magnitudeSquared();
        strainRateSq[cellIdx] = S(2.0) * symmMagSq;
    }

    return strainRateSq;
}


ScalarField RANS::velocityDivergence
(
    const FaceFluxField& flowRateFace
) const
{
    const Count numCells = mesh_.numDomainCells();

    ScalarField divU{mesh_};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const auto& cell = mesh_.cells()[cellIdx];
        const auto& faceIndices = cell.faceIndices();
        const auto& faceSigns = cell.faceSigns();

        Scalar sum = S(0.0);
        for (Index j = 0; j < faceIndices.size(); ++j)
        {
            sum += S(faceSigns[j]) * flowRateFace[faceIndices[j]];
        }

        divU[cellIdx] = sum / cell.volume();
    }

    return divU;
}


void RANS::updateResiduals
(
    const ScalarField& dissipationField,
    const ScalarField& dissipationPrevField
)
{
    lastKResidual_ = normalisedFieldResidual(k_, kPrev_);
    lastDissipationResidual_ =
        normalisedFieldResidual(dissipationField, dissipationPrevField);
}


Scalar RANS::normalisedFieldResidual
(
    const ScalarField& field,
    const ScalarField& previousField
) const
{
    const Count numCells = mesh_.numDomainCells();

    // Normalised change: ||x - x_prev|| / ||x_prev||
    Scalar diffSq = S(0.0);
    Scalar prevSq = S(0.0);

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar dField = field[cellIdx] - previousField[cellIdx];
        diffSq += dField * dField;
        prevSq += previousField[cellIdx] * previousField[cellIdx];
    }

    diffSq = globalSum(diffSq);
    prevSq = globalSum(prevSq);

    return
        std::sqrt(diffSq + vSmallValue)
      / std::sqrt(prevSq + vSmallValue);
}
