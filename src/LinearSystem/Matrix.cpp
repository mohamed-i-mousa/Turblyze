/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Matrix.cpp
 * @brief PETSc COO matrix assembly for finite volume transport equations
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Matrix.h"

// Standard library headers
#include <algorithm>
#include <cmath>
#include <vector>

// Project headers
#include "Comm.h"
#include "ErrorHandler.h"
#include "GlobalIndex.h"
#include "PETScRuntime.h"
#include "TimeScheme.h"
#include "HaloExchange.h"

// ************************* Special Member Functions *************************

Matrix::Matrix
(
    const Mesh& mesh,
    const BoundaryConditions& boundaryConds
)
:
    mesh_{mesh},
    bcManager_{boundaryConds}
{
    const Count numDomainCells = mesh_.numDomainCells();
    const Count numFaces = mesh_.numFaces();

    // Slot-ordered face lists; a ghost cell on either side means a cut
    faceSlot_.assign(numFaces, 0);

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            // A boundary face carries no COO slot: it hits only a diagonal
            boundaryFaces_.push_back(faceIdx);
        }
        else if
        (
            face.ownerCell() >= numDomainCells
         || face.neighborCell().value() >= numDomainCells
        )
        {
            processorFaces_.push_back(faceIdx);
        }
        else
        {
            faceSlot_[faceIdx] = 2 * internalFaces_.size();
            internalFaces_.push_back(faceIdx);
        }
    }

    const Count numInternalFaces = internalFaces_.size();
    const Count numProcessorFaces = processorFaces_.size();

    // COO layout: [internal pairs | processor | diagonals | scribble]
    for (Index m = 0; m < numProcessorFaces; ++m)
    {
        faceSlot_[processorFaces_[m]] = 2 * numInternalFaces + m;
    }

    diagOffset_ = 2 * numInternalFaces + numProcessorFaces;
    const Count numCoo = diagOffset_ + numDomainCells;
    scribbleSlot_ = numCoo;

    // The fixed sparsity pattern, in GLOBAL rank-major indices
    const GlobalIndex globalCells(numDomainCells);
    const IndexList& ghostGlobals = Halo::ghostGlobalIndices();

    std::vector<PetscInt> cooRows(numCoo);
    std::vector<PetscInt> cooCols(numCoo);

    for (Index k = 0; k < numInternalFaces; ++k)
    {
        const Face& face = mesh_.faces()[internalFaces_[k]];
        const auto owner =
            static_cast<PetscInt>(globalCells.toGlobal(face.ownerCell()));
        const auto neighbor =
            static_cast<PetscInt>
            (
                globalCells.toGlobal(face.neighborCell().value())
            );

        cooRows[2 * k] = owner;
        cooCols[2 * k] = neighbor;
        cooRows[2 * k + 1] = neighbor;
        cooCols[2 * k + 1] = owner;
    }

    // A processor face contributes one entry: owned row, ghost column
    for (Index m = 0; m < numProcessorFaces; ++m)
    {
        const Face& face = mesh_.faces()[processorFaces_[m]];
        const Index owner = face.ownerCell();
        const Index neighbor = face.neighborCell().value();
        const Index ownedSide = owner < numDomainCells ? owner : neighbor;
        const Index ghostSide = owner < numDomainCells ? neighbor : owner;

        cooRows[2 * numInternalFaces + m] =
            static_cast<PetscInt>(globalCells.toGlobal(ownedSide));
        cooCols[2 * numInternalFaces + m] =
            static_cast<PetscInt>
            (
                ghostGlobals[ghostSide - numDomainCells]
            );
    }

    for (Index cellIdx = 0; cellIdx < numDomainCells; ++cellIdx)
    {
        const auto diagIdx =
            static_cast<PetscInt>(globalCells.toGlobal(cellIdx));

        cooRows[diagOffset_ + cellIdx] = diagIdx;
        cooCols[diagOffset_ + cellIdx] = diagIdx;
    }

    const auto n = static_cast<PetscInt>(numDomainCells);

    CheckPETSc(MatCreate(PETScRuntime::comm(), &matrixA_));
    CheckPETSc(MatSetSizes(matrixA_, n, n, PETSC_DECIDE, PETSC_DECIDE));
    CheckPETSc(MatSetType(matrixA_, MATAIJ));
    CheckPETSc
    (
        MatSetPreallocationCOO
        (
            matrixA_,
            static_cast<PetscCount>(numCoo),
            cooRows.data(),
            cooCols.data()
        )
    );

    // COO value array with its processor-face scribble tail
    cooValues_.assign(numCoo + numProcessorFaces, S(0.0));
    vectorB_.assign(numDomainCells, S(0.0));

    // RHS handle wraps vectorB_ storage: staging writes are the vector
    CheckPETSc(VecCreate(PETScRuntime::comm(), &rhsVec_));
    CheckPETSc(VecSetSizes(rhsVec_, n, PETSC_DECIDE));
    CheckPETSc(VecSetType(rhsVec_, VECSTANDARD));
    CheckPETSc(VecPlaceArray(rhsVec_, vectorB_.data()));
}


Matrix::~Matrix() noexcept
{
    // Restore vectorB_ before destroying the PETSc vector
    if (rhsVec_ != nullptr)
    {
        VecResetArray(rhsVec_);
        VecDestroy(&rhsVec_);
    }

    if (matrixA_ != nullptr)
    {
        MatDestroy(&matrixA_);
    }
}

// ****************************** Public Methods ******************************

void Matrix::buildMatrix(const TransportEquation& equation)
{
    lastRelaxationFactor_ = S(0.0);

    const Count numDomainCells = mesh_.numDomainCells();
    const Count numInternalFaces = internalFaces_.size();
    const Count numProcessorFaces = processorFaces_.size();
    const Count numBoundaryFaces = boundaryFaces_.size();

    // Domain rows start from the cell source and any transient term; the
    // face loops accumulate on top, so this must run before all of them
    for (Index cellIdx = 0; cellIdx < numDomainCells; ++cellIdx)
    {
        Scalar diag = S(0.0);
        Scalar rhs = equation.source[cellIdx];

        // Transient term: add d(phi)/dt to the diagonal and RHS
        if (equation.transient)
        {
            const TransientTerm& transient = *equation.transient;

            const Scalar ddtPrevStep =
                transient.ddtPrevStep != nullptr
              ? (*transient.ddtPrevStep)[cellIdx]
              : S(0.0);

            const TimeContribution contribution =
                transient.scheme.coefficients
                (
                    mesh_.cells()[cellIdx].volume(),
                    transient.deltaT,
                    transient.phiPrevStep[cellIdx],
                    ddtPrevStep
                );

            diag += contribution.diag;
            rhs += contribution.source;
        }

        cooValues_[diagOffset_ + cellIdx] = diag;
        vectorB_[cellIdx] = rhs;
    }

    // Internal faces: each one owns its off-diagonal pair 2k / 2k + 1
    for (Index k = 0; k < numInternalFaces; ++k)
    {
        assembleInternalFace
        (
            2 * k,
            2 * k + 1,
            mesh_.faces()[internalFaces_[k]],
            equation
        );
    }

    // Same math across a cut, but only the domain cell's row has a slot
    for (Index m = 0; m < numProcessorFaces; ++m)
    {
        const Face& face = mesh_.faces()[processorFaces_[m]];
        const bool ownerIsLocal = face.ownerCell() < numDomainCells;
        const Index cooSlot = faceSlot_[processorFaces_[m]];

        assembleInternalFace
        (
            ownerIsLocal ? cooSlot : scribbleSlot_ + m,
            ownerIsLocal ? scribbleSlot_ + m : cooSlot,
            face,
            equation
        );
    }

    const BoundaryPatch* currentPatch = nullptr;
    const BoundaryType* currentBc = nullptr;

    for (Index m = 0; m < numBoundaryFaces; ++m)
    {
        const Face& face = mesh_.faces()[boundaryFaces_[m]];
        if (face.patch() != currentPatch)
        {
            currentPatch = face.patch();
            currentBc = &bcManager_.boundaryType(currentPatch->name(), equation.field);
        }
        assembleBoundaryFace(face, equation, *currentBc);
    }
}


void Matrix::relax(Scalar alpha, const ScalarField& phiPrevIter)
{
    if (alpha <= S(0.0))
    {
        FatalError("Matrix::relax: alpha must be positive");
    }

    const Count numCells = mesh_.numDomainCells();

    if (phiPrevIter.size() < numCells)
    {
        FatalError
        (
            "Matrix::relax: phiPrevIter size mismatch "
            "with matrix size"
        );
    }

    // Stored so setValues can recover the pre-relaxation diagonal
    lastRelaxationFactor_ = alpha;

    const Scalar factor = (S(1.0) - alpha) / alpha;

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar origDiag = cooValues_[diagOffset_ + cellIdx];

        // Scale diagonal: a_P <- a_P / alpha
        cooValues_[diagOffset_ + cellIdx] = origDiag / alpha;

        // Update RHS: b <- b + ((1-alpha)/alpha) * a_P * phiPrevIter
        vectorB_[cellIdx] += factor * origDiag * phiPrevIter[cellIdx];
    }
}


void Matrix::setValues
(
    const IndexList& cellIndices,
    const ScalarList& values,
    const ScalarField& ghostFractions,
    const ScalarField& ghostValues,
    const ScalarList& fractions
)
{
    const Count numDomainCells = mesh_.numDomainCells();
    const bool hasFractions = !fractions.empty();

    for (Index i = 0; i < cellIndices.size(); ++i)
    {
        const Index cellIdx = cellIndices[i];
        const Scalar f = hasFractions ? fractions[i] : S(1.0);

        if (f < rootSmallValue_)
        {
            continue;
        }

        const Scalar diag = cooValues_[diagOffset_ + cellIdx];

        if (f > S(1.0) - rootSmallValue_)
        {
            for (const Index faceIdx : mesh_.cells()[cellIdx].faceIndices())
            {
                const Face& face = mesh_.faces()[faceIdx];

                if (face.isBoundary())
                {
                    continue;
                }

                const Index slot = faceSlot_[faceIdx];
                const bool isOwner = face.ownerCell() == cellIdx;
                const Index neighborIdx =
                    isOwner ? face.neighborCell().value() : face.ownerCell();

                if (neighborIdx >= numDomainCells)
                {
                    // Processor face: the neighbor rank moves its own side
                    cooValues_[slot] = S(0.0);
                    continue;
                }

                const Index rowSlot = isOwner ? slot : slot + 1;
                const Index colSlot = isOwner ? slot + 1 : slot;

                const Scalar coupling = cooValues_[colSlot];
                if (coupling != S(0.0))
                {
                    vectorB_[neighborIdx] -= coupling * values[i];
                    cooValues_[colSlot] = S(0.0);
                }

                cooValues_[rowSlot] = S(0.0);
            }

            vectorB_[cellIdx] = diag * values[i];
        }
        else
        {
            const Scalar diagPre =
                lastRelaxationFactor_ > S(0.0)
              ? lastRelaxationFactor_ * diag
              : diag;

            const Scalar coeff = f / (S(1.0) - f) * diagPre;

            cooValues_[diagOffset_ + cellIdx] += coeff;
            vectorB_[cellIdx] += coeff * values[i];
        }
    }

    // Ghosts constrained elsewhere: move this rank's coupling to the RHS
    for (const Index faceIdx : processorFaces_)
    {
        const Face& face = mesh_.faces()[faceIdx];
        const Index owner = face.ownerCell();
        const Index neighbor = face.neighborCell().value();
        const Index ownedSide = owner < numDomainCells ? owner : neighbor;
        const Index ghostSide = owner < numDomainCells ? neighbor : owner;

        if (ghostFractions[ghostSide] <= S(1.0) - rootSmallValue_)
        {
            continue;
        }

        const Index slot = faceSlot_[faceIdx];
        const Scalar coupling = cooValues_[slot];

        if (coupling != S(0.0))
        {
            vectorB_[ownedSide] -= coupling * ghostValues[ghostSide];
            cooValues_[slot] = S(0.0);
        }
    }
}


void Matrix::assemble()
{
    CheckPETSc(MatSetValuesCOO(matrixA_, cooValues_.data(), INSERT_VALUES));

    PetscScalar* rhsArray = nullptr;
    CheckPETSc(VecGetArray(rhsVec_, &rhsArray));
    CheckPETSc(VecRestoreArray(rhsVec_, &rhsArray));
}


void Matrix::explicitJacobiUpdate
(
    const ScalarField& phiOld,
    ScalarField& phiNew
) const
{
    const Count numDomainCells = mesh_.numDomainCells();

    // phiNew_P = (b_P - sum_{N != P} A_PN phiOld_N) / A_PP
    for (Index cellIdx = 0; cellIdx < numDomainCells; ++cellIdx)
    {
        const Scalar diagonal = cooValues_[diagOffset_ + cellIdx];
        Scalar offDiagonalSum = S(0.0);

        for (const Index faceIdx : mesh_.cells()[cellIdx].faceIndices())
        {
            const Face& face = mesh_.faces()[faceIdx];

            if (face.isBoundary())
            {
                continue;
            }

            // Row cellIdx couples to the cell across each two-cell face
            const bool isOwner = face.ownerCell() == cellIdx;
            const Index neighborIdx =
                isOwner ? face.neighborCell().value() : face.ownerCell();
            const Index slot = faceSlot_[faceIdx];
            const Index rowSlot =
                neighborIdx >= numDomainCells
              ? slot
              : (isOwner ? slot : slot + 1);

            offDiagonalSum += cooValues_[rowSlot] * phiOld[neighborIdx];
        }

        phiNew[cellIdx] =
            (vectorB_[cellIdx] - offDiagonalSum) / (diagonal + vSmallValue);
    }
}

// ****************************** Private Methods *****************************

void Matrix::assembleInternalFace
(
    Index cooOwnerSlot,
    Index cooNeighborSlot,
    const Face& face,
    const TransportEquation& equation
)
{
    const Count numDomainCells = mesh_.numDomainCells();
    const Index ownerIdx = face.ownerCell();
    const Index neighborIdx = face.neighborCell().value();
    const Vector Sf = face.normal() * face.projectedArea();
    const Vector dPN = mesh_.dPN(face);
    const Scalar dPNMag = magnitude(dPN);
    const Vector ePN = dPN / (dPNMag + vSmallValue);

    // Orthogonal component (over-relaxed)
    const Vector Ef = (dot(Sf, Sf) / dot(Sf, ePN)) * ePN;
    const Scalar Gammaf = equation.GammaFace[face.idx()];
    const Scalar aDiff = Gammaf * magnitude(Ef) / (dPNMag + vSmallValue);

    // Convection coefficients
    Scalar flowRate = S(0.0);
    Scalar aPConv = S(0.0);
    Scalar aNConv = S(0.0);

    if (equation.convection)
    {
        const ConvectionTerm& convection = *equation.convection;

        flowRate = convection.flowRate[face.idx()];

        // Implicit upwind split: a_P = max(F, 0), a_N = min(F, 0)
        aPConv = std::max(flowRate, S(0.0));
        aNConv = std::min(flowRate, S(0.0));
    }

    // Off-diagonals own their slots, the diagonal shares scatter below
    cooValues_[cooOwnerSlot] = -aDiff + aNConv;     // A(owner, neighbor)
    cooValues_[cooNeighborSlot] = -aDiff - aPConv;  // A(neighbor, owner)

    // Non-orthogonal correction (explicit)
    const Vector Tf = Sf - Ef;
    const Vector& gradPhiP = equation.gradPhi[ownerIdx];
    const Vector& gradPhiN = equation.gradPhi[neighborIdx];
    const Vector gradPhif =
        equation.gradScheme.faceGradient
        (
            equation.phi,
            gradPhiP,
            gradPhiN,
            face.idx()
        );
    const Scalar nonOrthogonalFlux = Gammaf * dot(gradPhif, Tf);

    Scalar rhsOwner = nonOrthogonalFlux;
    Scalar rhsNeighbor = -nonOrthogonalFlux;

    // Deferred correction (explicit, convection only)
    if (equation.convection)
    {
        const ConvectionTerm& convection = *equation.convection;

        const Scalar deferredCorrection =
            convection.scheme.correction
            (
                mesh_,
                face,
                equation.phi,
                gradPhiP,
                gradPhiN,
                flowRate
            );

        rhsOwner -= deferredCorrection;
        rhsNeighbor += deferredCorrection;
    }

    // A ghost row belongs to the neighbor rank, so it takes nothing here
    if (ownerIdx < numDomainCells)
    {
        cooValues_[diagOffset_ + ownerIdx] += aDiff + aPConv;
        vectorB_[ownerIdx] += rhsOwner;
    }

    if (neighborIdx < numDomainCells)
    {
        cooValues_[diagOffset_ + neighborIdx] += aDiff - aNConv;
        vectorB_[neighborIdx] += rhsNeighbor;
    }
}


void Matrix::assembleBoundaryFace
(
    const Face& face,
    const TransportEquation& equation,
    const BoundaryType& bc
)
{
    const Index ownerIdx = face.ownerCell();
    const Index faceIdx = face.idx();
    const Index localIdx = faceIdx - face.patch()->firstFaceIdx();

    const Scalar GammaSf =
        equation.GammaFace[faceIdx] * face.projectedArea();

    const Scalar flux =
        equation.convection
      ? equation.convection->flowRate[faceIdx]
      : S(0.0);

    // Convection: flux * (a * phi_P + b)
    cooValues_[diagOffset_ + ownerIdx] += flux * bc.a(localIdx);
    vectorB_[ownerIdx]                 -= flux * bc.b(localIdx);

    // Diffusion: -GammaSf * (c * phi_P + d)
    cooValues_[diagOffset_ + ownerIdx] -= GammaSf * bc.c(localIdx);
    vectorB_[ownerIdx]                 += GammaSf * bc.d(localIdx);
}

