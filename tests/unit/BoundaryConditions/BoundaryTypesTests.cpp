/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryTypesTests.cpp
 * @brief Unit tests for the BoundaryType hierarchy and BoundaryConditions coefficients
 *
 * @details Tests BoundaryType property descriptors and BoundaryConditions
 * linearization coefficients (a, b, c, d) and face value reconstruction.
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <algorithm>
#include <memory>
#include <stdexcept>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "BoundaryType.h"
#include "FixedValue.h"
#include "ZeroGradient.h"
#include "FixedGradient.h"
#include "NoSlip.h"
#include "Symmetry.h"
#include "BoundaryConditions.h"
#include "BoundaryCoeffs.h"
#include "BoundaryTypeFactory.h"
#include "MeshFixtures.h"
#include "Field.h"
#include "Vector.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

const Scalar testOwnerValue = S(1.5);

/// Locate boundary patch by name
[[nodiscard]] const BoundaryPatch& findPatch(const Mesh& mesh, const Name& name)
{
    for (const auto& p : mesh.patches())
    {
        if (p.name() == name)
        {
            return p;
        }
    }
    throw std::runtime_error("patch not found: " + name);
}

/// Whether a name list contains a given token
[[nodiscard]] bool contains(const NameList& names, const Name& target)
{
    return std::find(names.begin(), names.end(), target) != names.end();
}

} // namespace

// ******************************** FixedValue ********************************

TEST_CASE("FixedValue linearization", "[bc]")
{
    FixedValue bc(S(2.0));

    REQUIRE(bc.typeName() == "fixedValue");
    REQUIRE(bc.value() == S(2.0));
    REQUIRE(bc.fixesValue());
    REQUIRE(!bc.isWall());
    REQUIRE(!bc.isSymmetry());

    const TestMesh box(1, 1, 1);
    const auto& patch = findPatch(box.mesh(), BoxPatch::xMin);
    const Index faceIdx = patch.firstFaceIdx();
    const Face& face = box.mesh().faces()[faceIdx];
    const Scalar dn = dot(box.mesh().dPf(face), face.normal());
    const Scalar diffMetric = S(1.0) / dn;

    bc.updateCoeffs(box.mesh(), patch);
    const Index localIdx = 0;

    // a = 0, b = V, c = -diffMetric, d = V * diffMetric
    REQUIRE_THAT(bc.a(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));
    REQUIRE_THAT(bc.b(localIdx), WithinRel(S(2.0), TestTolerances::relTight));
    REQUIRE_THAT(bc.c(localIdx), WithinRel(-diffMetric, TestTolerances::relTight));
    REQUIRE_THAT(bc.d(localIdx), WithinRel(S(2.0) * diffMetric, TestTolerances::relTight));

    // Reconstructed face value: a * phi_P + b = 2.0
    REQUIRE_THAT(bc.faceValue(localIdx, testOwnerValue), WithinRel(S(2.0), TestTolerances::relTight));

    BoundaryConditions::BCs bcs;
    auto bcPtr = std::make_unique<FixedValue>(S(2.0));
    bcPtr->updateCoeffs(box.mesh(), patch);
    bcs[BoxPatch::xMin][Field::p] = std::move(bcPtr);
    BoundaryConditions bcManager(std::move(bcs), box.mesh());

    const Scalar faceVal = bcManager.faceValue(face, testOwnerValue, Field::p);
    REQUIRE_THAT(faceVal, WithinRel(S(2.0), TestTolerances::relTight));
}

// ******************************* ZeroGradient *******************************

TEST_CASE("ZeroGradient linearization", "[bc]")
{
    ZeroGradient bc;

    REQUIRE(bc.typeName() == "zeroGradient");
    REQUIRE(!bc.fixesValue());
    REQUIRE(!bc.isWall());

    const TestMesh box(1, 1, 1);
    const auto& patch = findPatch(box.mesh(), BoxPatch::xMin);
    const Index faceIdx = patch.firstFaceIdx();
    const Face& face = box.mesh().faces()[faceIdx];

    bc.updateCoeffs(box.mesh(), patch);
    const Index localIdx = 0;

    // a = 1, b = 0, c = 0, d = 0
    REQUIRE_THAT(bc.a(localIdx), WithinRel(S(1.0), TestTolerances::relTight));
    REQUIRE_THAT(bc.b(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));
    REQUIRE_THAT(bc.c(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));
    REQUIRE_THAT(bc.d(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));

    // Reconstructed face value: a * phi_P + b = phi_P
    REQUIRE_THAT(bc.faceValue(localIdx, testOwnerValue), WithinRel(testOwnerValue, TestTolerances::relTight));

    BoundaryConditions::BCs bcs;
    auto bcPtr = std::make_unique<ZeroGradient>();
    bcPtr->updateCoeffs(box.mesh(), patch);
    bcs[BoxPatch::xMin][Field::p] = std::move(bcPtr);
    BoundaryConditions bcManager(std::move(bcs), box.mesh());

    const Scalar faceVal = bcManager.faceValue(face, testOwnerValue, Field::p);
    REQUIRE_THAT(faceVal, WithinRel(testOwnerValue, TestTolerances::relTight));
}

// ******************************* FixedGradient ******************************

TEST_CASE("FixedGradient linearization", "[bc]")
{
    FixedGradient bc(S(0.5));

    REQUIRE(bc.typeName() == "fixedGradient");
    REQUIRE(bc.gradient() == S(0.5));
    REQUIRE(bc.correctsBoundaryFlux());
    REQUIRE(!bc.fixesValue());

    const TestMesh box(1, 1, 1);
    const auto& patch = findPatch(box.mesh(), BoxPatch::xMin);
    const Index faceIdx = patch.firstFaceIdx();
    const Face& face = box.mesh().faces()[faceIdx];
    const Scalar dn = dot(box.mesh().dPf(face), face.normal());

    bc.updateCoeffs(box.mesh(), patch);
    const Index localIdx = 0;

    // a = 1, b = g * dn, c = 0, d = g
    REQUIRE_THAT(bc.a(localIdx), WithinRel(S(1.0), TestTolerances::relTight));
    REQUIRE_THAT(bc.b(localIdx), WithinRel(S(0.5) * dn, TestTolerances::relTight));
    REQUIRE_THAT(bc.c(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));
    REQUIRE_THAT(bc.d(localIdx), WithinRel(S(0.5), TestTolerances::relTight));

    // Reconstructed face value: a * phi_P + b = phi_P + g * dn
    REQUIRE_THAT(bc.faceValue(localIdx, testOwnerValue), WithinRel(testOwnerValue + S(0.5) * dn, TestTolerances::relTight));

    BoundaryConditions::BCs bcs;
    auto bcPtr = std::make_unique<FixedGradient>(S(0.5));
    bcPtr->updateCoeffs(box.mesh(), patch);
    bcs[BoxPatch::xMin][Field::p] = std::move(bcPtr);
    BoundaryConditions bcManager(std::move(bcs), box.mesh());

    const Scalar faceVal = bcManager.faceValue(face, testOwnerValue, Field::p);
    REQUIRE_THAT(faceVal, WithinRel(testOwnerValue + S(0.5) * dn, TestTolerances::relTight));
}

// ********************************** NoSlip **********************************

TEST_CASE("NoSlip is a zero-valued Dirichlet", "[bc]")
{
    NoSlip bc;

    REQUIRE(bc.typeName() == "noSlip");
    REQUIRE(bc.value() == S(0.0));
    REQUIRE(bc.fixesValue());
    REQUIRE(bc.isWall());

    const TestMesh box(1, 1, 1);
    const auto& patch = findPatch(box.mesh(), BoxPatch::xMin);
    const Index faceIdx = patch.firstFaceIdx();
    const Face& face = box.mesh().faces()[faceIdx];

    bc.updateCoeffs(box.mesh(), patch);
    const Index localIdx = 0;

    REQUIRE_THAT(bc.a(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));
    REQUIRE_THAT(bc.b(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));

    BoundaryConditions::BCs bcs;
    auto bcPtr = std::make_unique<NoSlip>();
    bcPtr->updateCoeffs(box.mesh(), patch);
    bcs[BoxPatch::xMin][Field::Ux] = std::move(bcPtr);
    BoundaryConditions bcManager(std::move(bcs), box.mesh());

    const Scalar faceVal = bcManager.faceValue(face, testOwnerValue, Field::Ux);
    REQUIRE_THAT(faceVal, WithinAbs(S(0.0), TestTolerances::absTight));
}

// ********************************* Symmetry *********************************

TEST_CASE("Symmetry linearization", "[bc]")
{
    Symmetry bcUx(Field::Ux);
    Symmetry bcUy(Field::Uy);

    REQUIRE(bcUx.typeName() == "symmetry");
    REQUIRE(bcUx.isSymmetry());

    const TestMesh box(1, 1, 1);
    const auto& patch = findPatch(box.mesh(), BoxPatch::xMin);
    const Index localIdx = 0;

    // Vector on xMin patch (normal = (-1, 0, 0), nx = -1, ny = 0, nz = 0)
    // a = 1 - nx^2 = 1 - (-1)^2 = 0
    ScalarField Ux(box.mesh(), S(2.0));
    ScalarField Uy(box.mesh(), S(3.0));
    ScalarField Uz(box.mesh(), S(0.0));

    bcUx.updateCoeffs(box.mesh(), patch);
    bcUx.refreshCoeffs(box.mesh(), patch, Ux, Uy, Uz);
    REQUIRE_THAT(bcUx.a(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));
    REQUIRE_THAT(bcUx.b(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));

    // For Uy on xMin (ny = 0, cross with nx * Ux):
    // a = 1 - ny^2 = 1 - 0 = 1
    // UnCross = nx * Ux + nz * Uz = -1 * 2.0 = -2.0
    // b = -ny * UnCross = 0
    bcUy.updateCoeffs(box.mesh(), patch);
    bcUy.refreshCoeffs(box.mesh(), patch, Ux, Uy, Uz);
    REQUIRE_THAT(bcUy.a(localIdx), WithinRel(S(1.0), TestTolerances::relTight));
    REQUIRE_THAT(bcUy.b(localIdx), WithinAbs(S(0.0), TestTolerances::absTight));
}

// ************************ BoundaryConditions faceValue **********************

TEST_CASE("BoundaryConditions::faceValue evaluates linear boundary reconstruction", "[bc]")
{
    TestMesh box(1, 1, 1);
    const auto& xMinPatch = findPatch(box.mesh(), BoxPatch::xMin);
    const auto& xMaxPatch = findPatch(box.mesh(), BoxPatch::xMax);

    BoundaryConditions::BCs bcs;
    auto fv = std::make_unique<FixedValue>(S(10.0));
    fv->updateCoeffs(box.mesh(), xMinPatch);
    auto zg = std::make_unique<ZeroGradient>();
    zg->updateCoeffs(box.mesh(), xMaxPatch);

    bcs[BoxPatch::xMin][Field::p] = std::move(fv);
    bcs[BoxPatch::xMax][Field::p] = std::move(zg);
    BoundaryConditions bcManager(std::move(bcs), box.mesh());

    const Face& xMinFace = box.mesh().faces()[xMinPatch.firstFaceIdx()];
    const Face& xMaxFace = box.mesh().faces()[xMaxPatch.firstFaceIdx()];

    REQUIRE_THAT(bcManager.faceValue(xMinFace, S(5.0), Field::p), WithinRel(S(10.0), TestTolerances::relTight));
    REQUIRE_THAT(bcManager.faceValue(xMaxFace, S(5.0), Field::p), WithinRel(S(5.0), TestTolerances::relTight));
}

// *************************** Selectable Type Names **************************

TEST_CASE("availableTypes lists the case-file-selectable BCs", "[bc]")
{
    const NameList velocityTypes = BoundaryTypeFactory::availableTypes(Field::Ux);

    REQUIRE(contains(velocityTypes, "fixedValue"));
    REQUIRE(contains(velocityTypes, "fixedGradient"));
    REQUIRE(contains(velocityTypes, "noSlip"));
    REQUIRE(contains(velocityTypes, "symmetry"));
    REQUIRE(contains(velocityTypes, "zeroGradient"));

    REQUIRE(contains(BoundaryTypeFactory::availableTypes(Field::p), "symmetry"));

    // The pressure correction field is not user-selectable
    REQUIRE(BoundaryTypeFactory::availableTypes(Field::pCorr).empty());
}