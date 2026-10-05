/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LinearSolverTests.cpp
 * @brief Both Krylov solvers reproduce a known solution on a tiny SPD system
 *
 * @details A 1D pure-diffusion problem on an 8x2x2 box decomposed across the
 * ranks. The exact solution is the linear profile phi(x) = x, which
 * finite-volume diffusion reproduces at the cell centres. The system is
 * assembled once and solved with each Krylov solver, so the profile must come
 * out the same for either solver at any rank count.
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <memory>
#include <span>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "Matrix.h"
#include "TransportEquation.h"
#include "LinearSolvers.h"
#include "BoundaryConditions.h"
#include "FixedValue.h"
#include "ZeroGradient.h"
#include "LeastSquares.h"
#include "MeshFixtures.h"
#include "CellData.h"
#include "FaceData.h"
#include "Field.h"
#include "Comm.h"
#include "StringTypes.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinAbs;

// ***************************** Internal Helpers *****************************

namespace
{

/// Register the fixed-value / zero-gradient BCs of the 1D diffusion problem
[[nodiscard]] BoundaryConditions makeDiffusionBoundaries(const Mesh& mesh)
{
    BoundaryConditions::BCs bcs;

    bcs[BoxPatch::xMin][Field::p] =
        std::make_unique<FixedValue>(S(0.0));
    bcs[BoxPatch::xMax][Field::p] =
        std::make_unique<FixedValue>(S(8.0));

    for
    (
        const Name& lateral :
        {
            BoxPatch::yMin,
            BoxPatch::yMax,
            BoxPatch::zMin,
            BoxPatch::zMax
        }
    )
    {
        bcs[lateral][Field::p] =
            std::make_unique<ZeroGradient>();
    }

    return BoundaryConditions(std::move(bcs), mesh);
}

} // namespace

// *********************** Krylov Solvers On A Box **************************

TEST_CASE("Krylov solvers reproduce the linear profile", "[petsc]")
{
    DecomposedBoxMesh box(8, 2, 2);
    const BoundaryConditions bc = makeDiffusionBoundaries(box.mesh());

    const LeastSquares gradScheme(box.mesh(), bc);

    ScalarField phi(box.mesh());
    const FaceFluxField gammaFace(box.mesh(), S(1.0));
    const ScalarField source(box.mesh());
    const VectorField gradPhi(box.mesh());

    const TransportEquation equation
    {
        .field = Field::p,
        .phi = phi,
        .GammaFace = gammaFace,
        .source = source,
        .gradPhi = gradPhi,
        .gradScheme = gradScheme
    };

    Matrix matrix(box.mesh(), bc);
    matrix.buildMatrix(equation);
    matrix.assemble();

    const std::vector<std::pair<Name, Name>> configurations =
    {
        {"PCG", "Jacobi"},
        {"PCG", "AMG"},
        {"PCG", "None"},
        {"BiCGSTAB", "Jacobi"},
        {"BiCGSTAB", "ILU"},
        {"BiCGSTAB", "BlockJacobi"},
        {"BiCGSTAB", "SOR"},
        {"BiCGSTAB", "AMG"},
        {"GMRES", "Jacobi"},
        {"GMRES", "ILU"},
        {"GMRES", "BlockJacobi"},
        {"GMRES", "AMG"}
    };

    for (const auto& [solverName, pcName] : configurations)
    {
        // A fresh zero-initialised solution vector per solver
        ScalarField solution(box.mesh());
        std::span<Scalar> x(solution.data(), box.mesh().numDomainCells());

        const auto solver = LinearSolver::create
        (
            solverName,
            pcName,
            TestTolerances::solverTolerance,
            Count{200},
            "test"
        );
        solver->solve(x, matrix.matrixA(), matrix.rhsVec());

        INFO("Testing solver=" << solverName << ", preconditioner=" << pcName);
        CHECK(solver->lastPerformance().converged);
        CHECK(solver->lastPerformance().solverName == solverName);
        CHECK(solver->preconditioner() == pcName);

        for
        (
            Index cellIdx = 0;
            cellIdx < box.mesh().numDomainCells();
            ++cellIdx
        )
        {
            CHECK_THAT
            (
                solution[cellIdx],
                WithinAbs
                (
                    box.mesh().cells()[cellIdx].centroid().x(),
                    TestTolerances::absSolve
                )
            );
        }
    }
}