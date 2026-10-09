/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PotentialFlow.cpp
 * @brief Implementation of potential flow initialization via Poisson eq.
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "PotentialFlow.h"

// External library headers
#include <petscmat.h>

// Project headers
#include "Matrix.h"
#include "LinearSolvers.h"
#include "PETScRuntime.h"
#include "ErrorHandler.h"
#include "HaloExchange.h"
#include "Reduce.h"
#include "TransportEquation.h"
#include "BoundaryPatch.h"
#include "Face.h"
#include "GradientScheme.h"

// ************************* namespace Initialization *************************

namespace Initialization
{

// *************************** Initialization Method **************************

void PotentialFlow::initialize
(
    const Mesh& mesh,
    BoundaryConditions& bc,
    ScalarField& Ux,
    ScalarField& Uy,
    ScalarField& Uz,
    ScalarField& p
) const
{
    const Count numCells = mesh.numDomainCells();

    // Initialize velocity with freestream
    Ux.setAll(Uinf_.x());
    Uy.setAll(Uinf_.y());
    Uz.setAll(Uinf_.z());

    // Initial face fluxes
    FaceFluxField flowRateFace(mesh, S(0.0));

    for (const Face& f : mesh.faces())
    {
        if (!f.isBoundary())
        {
            flowRateFace[f.idx()] =
                dot(Uinf_, f.normal() * f.projectedArea());
        }
    }

    for (const BoundaryPatch& patch : mesh.patches())
    {
        if (patch.type() == PatchType::processor)
        {
            continue;
        }

        for
        (
            Index faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            const Face& f = mesh.faces()[faceIdx];
            const Index owner = f.ownerCell();
            const Vector Uf
            (
                bc.faceValue(f, Ux[owner], Field::Ux),
                bc.faceValue(f, Uy[owner], Field::Uy),
                bc.faceValue(f, Uz[owner], Field::Uz)
            );
            flowRateFace[faceIdx] =
                bc.boundaryType(f, Field::Ux).isSymmetry()
              ? S(0.0)
              : dot(Uf, f.normal() * f.projectedArea());
        }
    }

    // Set up fields
    ScalarField massImbalance(mesh, S(0.0));
    ScalarField phi(mesh, S(0.0));
    VectorField gradPhi(mesh, Vector{S(0.0), S(0.0), S(0.0)});
    FaceFluxField GammaFace(mesh, S(1.0));
    Matrix matrix(mesh, bc);

    // Solve velocity potential with Non-orthogonal corrector loops
    for (Count corr = 0; corr < numNonOrthoCorr_; ++corr)
    {
        for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            Scalar net = S(0.0);
            const auto& faceIndices = mesh.cells()[cellIdx].faceIndices();
            const auto& signs = mesh.cells()[cellIdx].faceSigns();

            for (Index j = 0; j < faceIndices.size(); ++j)
            {
                net += signs[j] * flowRateFace[faceIndices[j]];
            }

            massImbalance[cellIdx] = -net;
        }

        phi.setAll(S(0.0));

        TransportEquation poissonEq
        {
            .field      = Field::pCorr,
            .phi        = phi,
            .convection = std::nullopt,
            .GammaFace  = GammaFace,
            .source     = massImbalance,
            .gradPhi    = gradPhi,
            .gradScheme = gradScheme_
        };

        matrix.buildMatrix(poissonEq);
        matrix.assemble();

        solver_.solve
        (
            {phi.data(), numCells},
            matrix.matrixA(),
            matrix.rhsVec()
        );

        Halo::exchange({&phi});
        gradScheme_.fieldGradient(Field::pCorr, phi, gradPhi);
        Halo::exchange({&gradPhi});

        // Correct velocity
        for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            Ux[cellIdx] -= gradPhi[cellIdx].x();
            Uy[cellIdx] -= gradPhi[cellIdx].y();
            Uz[cellIdx] -= gradPhi[cellIdx].z();
        }
        Halo::exchange({&Ux, &Uy, &Uz});

        // Correct internal face fluxes
        for (const Face& f : mesh.faces())
        {
            if (!f.isBoundary())
            {
                const Index owner = f.ownerCell();
                const Index neighbor = f.neighborCell().value();
                const Vector gf =
                    gradScheme_.faceGradient
                    (
                        phi,
                        gradPhi[owner],
                        gradPhi[neighbor],
                        f.idx()
                    );
                flowRateFace[f.idx()] -=
                    dot(gf, f.normal() * f.projectedArea());
            }
        }
    }

    // Compute Bernoulli pressure: p = pinf + 0.5 * rho * (Uinf^2 - |U|^2)
    const Scalar UinfSq = magnitudeSquared(Uinf_);
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar USq =
            Ux[cellIdx] * Ux[cellIdx]
          + Uy[cellIdx] * Uy[cellIdx]
          + Uz[cellIdx] * Uz[cellIdx];

        p[cellIdx] = pinf_ + S(0.5) * rho_ * (UinfSq - USq);
    }

    Halo::exchange({&p});
}

} // namespace Initialization
