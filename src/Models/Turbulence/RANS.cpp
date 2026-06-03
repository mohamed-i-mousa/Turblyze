/******************************************************************************
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
#include "BoundaryData.h"
#include "BoundaryConditions.h"
#include "Field.h"
#include "Matrix.h"

// ************************* Special Member Functions *************************

RANS::RANS
(
    const Mesh& mesh,
    const BoundaryConditions& bc,
    const GradientScheme& gradScheme,
    const ConvectionSchemes& kScheme,
    LinearSolver& kSolver,
    const ConvectionSchemes& dissipationScheme,
    LinearSolver& dissipationSolver,
    const Scalar nu,
    const Scalar alphaK,
    const Scalar alphaDissipation,
    const bool debug
)
:
    mesh_{mesh},
    bcManager_{bc},
    gradientScheme_{gradScheme},
    matrixConstruct_{std::make_unique<Matrix>(mesh_, bcManager_)},
    kConvectionScheme_{kScheme},
    kSolver_{kSolver},
    dissipationConvectionScheme_{dissipationScheme},
    dissipationSolver_{dissipationSolver},
    nu_{nu},
    alphaK_{alphaK},
    alphaDissipation_{alphaDissipation},
    debug_{debug}
{}

RANS::~RANS() noexcept = default;

// ************************ Inlet Condition Calculators ***********************

Scalar RANS::inletK
(
    const Vector& velocity,
    const Scalar turbulenceIntensity
) noexcept
{
    const Scalar uPrime = turbulenceIntensity * magnitude(velocity);
    return std::max(S(1.5) * uPrime * uPrime, smallValue);
}

// ***************************** Accessor Methods *****************************

Scalar RANS::boundaryTurbulentViscosity
(
    const Face& face,
    const BoundaryConditions& bcManager
) const
{
    if (face.isBoundary())
    {
        const BoundaryPatch& patch = face.patch()->get();
        const BoundaryData& bc =
            bcManager.fieldBC(patch.patchName(), Field::nut);

        if (bc.type() == BCType::nutWallFunction)
        {
            return nutWall_[face.idx()];
        }
    }

    return TurbulenceModel::boundaryTurbulentViscosity(face, bcManager);
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
        {"yPlus", &yPlus_},
        {"wallShearStress", &wallShearStress_}
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

void RANS::updateWallDistance()
{
    wallDistanceConverged_ = false;
    wallDistance_.setAll(S(1e10));
    nearestWallPoint_.setAll(Vector{S(1e15), S(1e15), S(1e15)});

    // Seed wall-adjacent cells with the perpendicular distance to each
    // wall face centroid
    for (const auto& face : mesh_.faces())
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch& patch = face.patch()->get();

        if (patch.type() != PatchType::wall) continue;

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

        if (maxChange < tolerance)
        {
            wallDistanceConverged_ = true;
            break;
        }
    }
}


void RANS::updateYPlusLam(const Scalar kappa, const Scalar E)
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
    const Field wallFunctionField,
    const BCType wallFunctionType
)
{
    const Count numCells = mesh_.numCells();

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

        const BoundaryPatch& patch = face.patch()->get();
        if (patch.type() != PatchType::wall) continue;

        const BoundaryData& bc =
            bcManager.fieldBC(patch.patchName(), wallFunctionField);
        if (bc.type() != wallFunctionType) continue;

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
            std::abs(dot(face.dPf(), face.normal())),
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

            const auto& patch = face.patch();
            if
            (
                patch.has_value()
             && patch->get().type() == PatchType::wall
            )
            {
                totalPolyWallArea[cellIdx] += face.projectedArea();
            }
        }
    }

    wallCellFraction_.resize(wallCellIndices_.size());

    constexpr Scalar wallCellFractionTol = S(0.1);

    #pragma omp parallel for schedule(static)
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


void RANS::updateYPlus
(
    const ScalarField& turbulentKineticEnergy,
    const Scalar cmu25
)
{
    for (Index faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const Index cellIdx = face.ownerCell();

        yPlus_[face.idx()] =
            cmu25
          * std::sqrt(turbulentKineticEnergy[cellIdx])
          * y_[face.idx()]
          / nu_;
    }
}


void RANS::updateWallShearStress
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz,
    const ScalarField& turbulentKineticEnergy,
    const Scalar cmu25
)
{
    for (Index faceIdx : wallFunctionFaceIndices_)
    {
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
          : cmu25 * cmu25 * turbulentKineticEnergy[cellIdx];

        wallShearStress_[face.idx()] = std::min(tau, S(1000.0));
    }
}


ScalarField RANS::computeStrainRateMagnitude(const TensorField& gradU) const
{
    const Count numCells = mesh_.numCells();

    ScalarField strainRateMag;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // S = sqrt(2 * S_ij * S_ij), S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        const Scalar symmMagSq = gradU[cellIdx].symm().magnitudeSquared();
        strainRateMag[cellIdx] = std::sqrt(S(2.0) * symmMagSq);
    }

    return strainRateMag;
}


ScalarField RANS::velocityDivergence
(
    const FaceFluxField& flowRateFace
) const
{
    const Count numCells = mesh_.numCells();

    ScalarField divU;

    #pragma omp parallel for schedule(static)
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
    const ScalarField& dissipation,
    const ScalarField& dissipationPrev
)
{
    lastKResidual_ = normalisedFieldResidual(k_, kPrev_);
    lastDissipationResidual_ =
        normalisedFieldResidual(dissipation, dissipationPrev);
}


Scalar RANS::normalisedFieldResidual
(
    const ScalarField& field,
    const ScalarField& previousField
) const
{
    const Count numCells = mesh_.numCells();

    // Normalised change: ||x - x_prev|| / ||x_prev||
    Scalar diffSq = S(0.0);
    Scalar prevSq = S(0.0);

    #pragma omp parallel for schedule(static) reduction(+:diffSq, prevSq)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar dField = field[cellIdx] - previousField[cellIdx];
        diffSq += dField * dField;
        prevSq += previousField[cellIdx] * previousField[cellIdx];
    }

    return
        std::sqrt(diffSq + vSmallValue)
      / std::sqrt(prevSq + vSmallValue);
}
