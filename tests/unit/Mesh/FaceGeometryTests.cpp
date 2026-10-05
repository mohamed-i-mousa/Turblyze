/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FaceGeometryTests.cpp
 * @brief Unit tests for Face area, centroid, normal, and cell distances
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "MeshFixtures.h"
#include "TestTolerances.h"
#include "Mesh.h"
#include "Face.h"
#include "Vector.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// Return the first internal face, or nullptr if none exist
[[nodiscard]] const Face* firstInternalFace(const Mesh& mesh)
{
    for (Index faceIdx = 0; faceIdx < mesh.numFaces(); ++faceIdx)
    {
        const Face& face = mesh.faces()[faceIdx];

        if (!face.isBoundary())
        {
            return &face;
        }
    }

    return nullptr;
}

} // namespace

// **************************** Unit-Square Quad ******************************

TEST_CASE("Unit-square quad geometry", "[mesh]")
{
    // Four corners of a unit square in the z = 0 plane
    const NodeList nodes
    {
        Vector(S(0.0), S(0.0), S(0.0)),
        Vector(S(1.0), S(0.0), S(0.0)),
        Vector(S(1.0), S(1.0), S(0.0)),
        Vector(S(0.0), S(1.0), S(0.0))
    };

    Face face(0, IndexList{0, 1, 2, 3}, 0);
    static_cast<void>(face.geometricProperties(nodes));

    REQUIRE_THAT
    (
        face.projectedArea(),
        WithinRel(S(1.0), TestTolerances::relTight)
    );

    REQUIRE_THAT
    (
        face.centroid().x(),
        WithinAbs(S(0.5), TestTolerances::absTight)
    );
    REQUIRE_THAT
    (
        face.centroid().y(),
        WithinAbs(S(0.5), TestTolerances::absTight)
    );
    REQUIRE_THAT
    (
        face.centroid().z(),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );

    REQUIRE_THAT
    (
        magnitude(face.normal()),
        WithinRel(S(1.0), TestTolerances::relTight)
    );
}

// ******************************* Triangle Area ******************************

TEST_CASE("Triangle area", "[mesh]")
{
    // A right triangle with legs of length one: area = 0.5
    const NodeList nodes
    {
        Vector(S(0.0), S(0.0), S(0.0)),
        Vector(S(1.0), S(0.0), S(0.0)),
        Vector(S(0.0), S(1.0), S(0.0))
    };

    Face face(0, IndexList{0, 1, 2}, 0);
    static_cast<void>(face.geometricProperties(nodes));

    REQUIRE_THAT
    (
        face.projectedArea(),
        WithinRel(S(0.5), TestTolerances::relTight)
    );
}

// ************************** Two-Cell Box Distances **************************

TEST_CASE("Two-cell box face distances", "[mesh]")
{
    // Two unit cells along x share one internal face at x = 1. The owner and
    // neighbour centroids sit half a cell either side, so both distances = 0.5
    const TestMesh box(2, 1, 1);

    const Face* internalFace = firstInternalFace(box.mesh());
    REQUIRE(internalFace != nullptr);

    REQUIRE_THAT
    (
        magnitude(box.mesh().dPf(*internalFace)),
        WithinRel(S(0.5), TestTolerances::relTight)
    );

    REQUIRE_THAT
    (
        magnitude(box.mesh().dNf(*internalFace)),
        WithinRel(S(0.5), TestTolerances::relTight)
    );

    REQUIRE_THAT
    (
        magnitude(internalFace->normal()),
        WithinRel(S(1.0), TestTolerances::relTight)
    );
}

// ************************** Boundary Face Diffusion *************************

TEST_CASE("Boundary face geometric diffusion factor", "[mesh]")
{
    // Single unit cube: owner centroid at (0.5, 0.5, 0.5), face centroid at
    // boundary (e.g. x = 0), |dPf| = 0.5, orthogonal so |Ef| == |Sf| = 1.0.
    // gDiff = |Ef| / (|Sf| * |dPf|) = 1.0 / (1.0 * 0.5) = 2.0
    const TestMesh box(1, 1, 1);
    const Face& boundaryFace = box.mesh().faces().front();

    REQUIRE_THAT
    (
        box.mesh().gDiff(boundaryFace),
        WithinRel(S(2.0), TestTolerances::relTight)
    );
}