/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file GradientLimiterTests.cpp
 * @brief Barth-Jespersen gradient limiter on the interior of a hex box
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <algorithm>
#include <memory>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "MeshFixtures.h"
#include "TestTolerances.h"
#include "LeastSquares.h"
#include "BoundaryConditions.h"
#include "ZeroGradient.h"
#include "CellData.h"
#include "Vector.h"
#include "Field.h"

using Catch::Matchers::WithinAbs;

// ***************************** Internal Helpers *****************************

namespace
{

/// Interior centre cell of a 3x3x3 box (all six faces internal)
constexpr Index centreCell = 13;

/// The linear field phi = 1.5 x + 2.5 y - 3.5 z + 4.2
[[nodiscard]] Scalar linearField(const Vector& x) noexcept
{
    return S(1.5) * x.x() + S(2.5) * x.y() - S(3.5) * x.z() + S(4.2);
}

/// Register a zero-gradient pressure BC on every patch of the box
[[nodiscard]] BoundaryConditions makeZeroGradient(const Mesh& mesh)
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
        bcs[name][Field::p] = std::make_unique<ZeroGradient>();
    }

    return BoundaryConditions(std::move(bcs), mesh);
}

} // namespace

// ************************ Linear Field Passes Through **********************

TEST_CASE("The limiter leaves a linear-field gradient intact", "[schemes]")
{
    TestMesh box(3, 3, 3);
    const BoundaryConditions bc = makeZeroGradient(box.mesh());

    const LeastSquares leastSquares(box.mesh(), bc);

    ScalarField phi(box.mesh());
    for (Index cellIdx = 0; cellIdx < box.mesh().numCells(); ++cellIdx)
    {
        phi[cellIdx] = linearField(box.mesh().cells()[cellIdx].centroid());
    }

    // fieldGradient computes each cell gradient then limits it in place
    VectorField gradPhi(box.mesh());
    leastSquares.fieldGradient(Field::p, phi, gradPhi);

    // The interior cell reconstruction never overshoots, so alpha stays 1
    REQUIRE_THAT
    (
        gradPhi[centreCell].x(),
        WithinAbs(S(1.5),
        TestTolerances::absOperator)
    );
    
    REQUIRE_THAT
    (
        gradPhi[centreCell].y(),
        WithinAbs(S(2.5),
        TestTolerances::absOperator)
    );
    
    REQUIRE_THAT
    (
        gradPhi[centreCell].z(),
        WithinAbs(S(-3.5),
        TestTolerances::absOperator)
    );
}

// ************************ Over-Steep Gradient Clamped **********************

TEST_CASE("The limiter clamps an over-steep gradient to the hull", "[schemes]")
{
    TestMesh box(3, 3, 3);
    const BoundaryConditions bc = makeZeroGradient(box.mesh());

    const LeastSquares leastSquares(box.mesh(), bc);

    ScalarField phi(box.mesh());
    for (Index cellIdx = 0; cellIdx < box.mesh().numCells(); ++cellIdx)
    {
        phi[cellIdx] = linearField(box.mesh().cells()[cellIdx].centroid());
    }

    VectorField gradPhi(box.mesh());
    leastSquares.fieldGradient(Field::p, phi, gradPhi);

    // Inject a wildly over-steep gradient into the interior cell only
    const Vector overSteep = S(50.0) * gradPhi[centreCell];
    gradPhi[centreCell] = overSteep;

    leastSquares.limitGradient(Field::p, phi, gradPhi);

    // The neighbour value hull the reconstruction must stay inside
    const Cell& cell = box.mesh().cells()[centreCell];
    Scalar phiMin = phi[centreCell];
    Scalar phiMax = phi[centreCell];
    for (const Index neighborIdx : cell.neighborCellIndices())
    {
        phiMin = std::min(phiMin, phi[neighborIdx]);
        phiMax = std::max(phiMax, phi[neighborIdx]);
    }

    // Every face reconstruction now lies within the hull
    for (const Index faceIdx : cell.faceIndices())
    {
        const Vector r =
            box.mesh().faces()[faceIdx].centroid() - cell.centroid();
        const Scalar phiFace = phi[centreCell] + dot(gradPhi[centreCell], r);

        REQUIRE(phiFace <= phiMax + TestTolerances::absOperator);
        REQUIRE(phiFace >= phiMin - TestTolerances::absOperator);
    }

    // The clamp strictly reduced the injected gradient (alpha < 1)
    REQUIRE(magnitude(gradPhi[centreCell]) < magnitude(overSteep));
}