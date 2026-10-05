/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TurbulenceBehaviorTests.cpp
 * @brief k-omega SST construction and one solve on a flat-wall box
 *
 * @details A standalone kOmegaSST is stood up on a unit hex box whose six
 * patches are all walls. Construction runs the wall-distance wave, whose
 * exact answer on a box is the shortest gap to one of the six wall planes.
 * One solve with a hand-built shear then exercises the private blending /
 * production / bound kernels through public behaviour. The box is decomposed
 * across the ranks, so the wave and the solve both cross a cut.
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <algorithm>
#include <cmath>
#include <memory>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "kOmegaSST.h"
#include "TimeScheme.h"
#include "GradientScheme.h"
#include "ConvectionScheme.h"
#include "LinearSolvers.h"
#include "BoundaryConditions.h"
#include "ZeroGradient.h"
#include "NoSlip.h"
#include "MeshFixtures.h"
#include "CellData.h"
#include "FaceData.h"
#include "Tensor.h"
#include "Field.h"
#include "Comm.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// Register wall BCs (NoSlip for Ux, ZeroGradient for k, omega, nut) on every wall patch.
[[nodiscard]] BoundaryConditions makeTurbulenceBoundaries(const Mesh& mesh)
{
    BoundaryConditions::BCs bcs;

    for
    (
        const Name& name :
        {
            BoxPatch::xMin,
            BoxPatch::xMax,
            BoxPatch::yMin,
            BoxPatch::yMax,
            BoxPatch::zMin,
            BoxPatch::zMax
        }
    )
    {
        for (const Field field : {Field::Ux, Field::Uy, Field::Uz})
        {
            bcs[name][field] = std::make_unique<NoSlip>();
        }

        for (const Field field : {Field::k, Field::omega, Field::nut})
        {
            bcs[name][field] = std::make_unique<ZeroGradient>();
        }
    }

    return BoundaryConditions(std::move(bcs), mesh);
}

/// Shortest distance from a point to the six faces of the cubic box
[[nodiscard]] Scalar planeWallDistance(const Vector& point, Scalar side)
{
    const Scalar toMin =
        std::min({point.x(), point.y(), point.z()});
    const Scalar toMax =
        std::min({side - point.x(), side - point.y(), side - point.z()});

    return std::min(toMin, toMax);
}

/// Locate a named cell-centred output field of a turbulence model
[[nodiscard]] const ScalarField* findOutput
(
    const TurbulenceModel::CellDataPair& outputs,
    const Name& name
)
{
    for (const auto& [outputName, field] : outputs)
    {
        if (outputName == name)
        {
            return field;
        }
    }
    return nullptr;
}

} // namespace

// ******************** Construction And One Solve **********************

TEST_CASE
(
    "kOmegaSST builds wall distance and bounds k/omega",
    "[petsc][models]"
)
{
    DecomposedBoxMesh box(3, 3, 3);
    const BoundaryConditions bc = makeTurbulenceBoundaries(box.mesh());

    // Standalone dependencies; every one must outlive the model
    const auto timeScheme = TimeScheme::create("steadyState");
    const auto gradScheme =
        GradientScheme::create("leastSquares", box.mesh(), bc);
    const auto kScheme = ConvectionScheme::create("Upwind");
    const auto omegaScheme = ConvectionScheme::create("Upwind");
    const auto kSolver = LinearSolver::create
    (
        "PCG",
        "Jacobi",
        TestTolerances::solverTolerance,
        Count{200},
        "k"
    );
    const auto omegaSolver = LinearSolver::create
    (
        "PCG",
        "Jacobi",
        TestTolerances::solverTolerance,
        Count{200},
        "omega"
    );

    const Scalar initialK = S(0.375);
    const Scalar initialOmega = S(100.0);

    kOmegaSST model
    (
        box.mesh(),
        bc,
        *timeScheme,
        *gradScheme,
        *kScheme,
        *kSolver,
        *omegaScheme,
        *omegaSolver,
        S(1.0),          // deltaT (steady: unused)
        S(1.0e-5),       // nu
        initialK,
        initialOmega,
        S(0.7),          // alphaK
        S(0.7),          // alphaOmega
        false,           // roughWall
        false            // debug
    );

    // A const view: k() has a protected mutable overload, so the reads below
    // must resolve to the public const one (solve() still uses the model)
    const kOmegaSST& view = model;

    // Wall distance on the flat-wall box. Everything checked ahead of the
    // collective solve() uses CHECK: REQUIRE would unwind the failing rank
    // out of the test case and leave the others waiting inside solve()

    CHECK(view.wallDistanceConverged());

    const ScalarField* wallDistance =
        findOutput(view.cellDataOutputs(), "wallDistance");
    CHECK(wallDistance != nullptr);

    // Every wall is a box face, so the exact distance is the shortest gap to
    // one of the six planes - checked per owned cell, at any rank count
    if (wallDistance != nullptr)
    {
        for (Index cellIdx = 0; cellIdx < box.mesh().numDomainCells(); ++cellIdx)
        {
            CHECK_THAT
            (
                (*wallDistance)[cellIdx],
                WithinAbs
                (
                    planeWallDistance
                    (
                        box.mesh().cells()[cellIdx].centroid(),
                        S(3.0)
                    ),
                    TestTolerances::absOperator
                )
            );
        }
    }

    // Initial field values

    for (Index cellIdx = 0; cellIdx < box.mesh().numDomainCells(); ++cellIdx)
    {
        CHECK_THAT
        (
            view.k()[cellIdx],
            WithinRel(initialK, TestTolerances::relTight)
        );
        CHECK_THAT
        (
            view.dissipation()[cellIdx],
            WithinRel(initialOmega, TestTolerances::relTight)
        );
        CHECK_THAT
        (
            view.turbulentViscosity()[cellIdx],
            WithinRel(initialK / initialOmega, TestTolerances::relTight)
        );
    }

    // One solve with a hand-built shear

    const ScalarField Ux(box.mesh());
    const ScalarField Uy(box.mesh());
    const ScalarField Uz(box.mesh());
    const FaceFluxField flowRateFace(box.mesh(), S(0.0));

    // Pure shear du_x/dy = 2, so the strain-rate magnitude is nonzero
    TensorField gradU(box.mesh());
    for (Index cellIdx = 0; cellIdx < box.mesh().numCells(); ++cellIdx)
    {
        gradU[cellIdx] = Tensor
        (
            S(0.0), S(2.0), S(0.0),
            S(0.0), S(0.0), S(0.0),
            S(0.0), S(0.0), S(0.0)
        );
    }

    model.solve(Ux, Uy, Uz, flowRateFace, gradU);

    // k and omega stay strictly positive; nut stays finite and non-negative
    for (Index cellIdx = 0; cellIdx < box.mesh().numDomainCells(); ++cellIdx)
    {
        REQUIRE(view.k()[cellIdx] > S(0.0));
        REQUIRE(view.dissipation()[cellIdx] > S(0.0));
        REQUIRE(std::isfinite(view.turbulentViscosity()[cellIdx]));
        REQUIRE(view.turbulentViscosity()[cellIdx] >= S(0.0));
    }

    // The residuals are updated to finite, non-negative values
    for (const auto& [name, value] : view.residualOutputs())
    {
        REQUIRE(std::isfinite(value));
        REQUIRE(value >= S(0.0));
    }
}