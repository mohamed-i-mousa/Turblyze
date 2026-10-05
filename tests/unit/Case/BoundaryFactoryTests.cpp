/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryFactoryTests.cpp
 * @brief BoundaryType::create from case sections
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <memory>
#include <string>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "BoundaryType.h"
#include "BoundaryTypeFactory.h"
#include "FixedValue.h"
#include "FixedGradient.h"
#include "ZeroGradient.h"
#include "NoSlip.h"
#include "Symmetry.h"
#include "CaseReader.h"
#include "Field.h"
#include "Vector.h"
#include "StringTypes.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{



/// The boundaryTypes section of the committed parser fixture
[[nodiscard]] CaseReader boundaryTypesFixture()
{
    const CaseReader reader
    (
        FilePath(TURBLYZE_TEST_FIXTURE_DIR) + "/cases/parserCase"
    );
    return reader.section("boundaryTypes");
}

} // namespace

// ************************* Fixed-Value Velocity ****************************

TEST_CASE("Factory builds a fixedValue velocity component", "[bc][selection]")
{
    const CaseReader bcSection = boundaryTypesFixture();
    const CaseReader& patch = bcSection.section("fixedVelocity");
    const Vector val = patch.lookup<Vector>("value");

    // Each velocity component draws its own scalar from the (1 2 3) vector
    const auto bcX =
        BoundaryTypeFactory::create("fixedValue", Field::Ux, val.x());
    const auto bcZ =
        BoundaryTypeFactory::create("fixedValue", Field::Uz, val.z());

    REQUIRE(bcX->typeName() == "fixedValue");
    REQUIRE(bcX->fixesValue());
    REQUIRE(dynamic_cast<const FixedValue*>(bcX.get())->value() == S(1.0));
    REQUIRE(dynamic_cast<const FixedValue*>(bcZ.get())->value() == S(3.0));
}

// ************************** Fixed-Value Scalar *****************************

TEST_CASE("Factory builds a fixedValue scalar", "[bc][selection]")
{
    const CaseReader bcSection = boundaryTypesFixture();
    const CaseReader& patch = bcSection.section("fixedScalar");
    const Scalar val = patch.lookup<Scalar>("value");

    const auto bc = BoundaryTypeFactory::create
    (
        "fixedValue", Field::p, val
    );

    REQUIRE(bc->typeName() == "fixedValue");
    REQUIRE(dynamic_cast<const FixedValue*>(bc.get())->value() == S(2.5));
}

// ************************** Fixed-Gradient Scalar **************************

TEST_CASE("Factory builds a fixedGradient scalar", "[bc][selection]")
{
    const CaseReader bcSection = boundaryTypesFixture();
    const CaseReader& patch = bcSection.section("gradientScalar");
    const Scalar grad = patch.lookup<Scalar>("gradient");

    const auto bc = BoundaryTypeFactory::create
    (
        "fixedGradient", Field::p, grad
    );

    REQUIRE(bc->typeName() == "fixedGradient");
    REQUIRE(bc->correctsBoundaryFlux());
    REQUIRE(dynamic_cast<const FixedGradient*>(bc.get())->gradient() == S(0.5));
}

// ****************************** Trait Types ********************************

TEST_CASE("Factory builds the zero-parameter trait types", "[bc][selection]")
{
    const auto slip = BoundaryTypeFactory::create
    (
        "zeroGradient", Field::p
    );
    REQUIRE(slip->typeName() == "zeroGradient");

    const auto wall = BoundaryTypeFactory::create
    (
        "noSlip", Field::Ux
    );
    REQUIRE(wall->typeName() == "noSlip");
    REQUIRE(wall->fixesValue());
    REQUIRE(dynamic_cast<const NoSlip*>(wall.get())->value() == S(0.0));
}

// ******************************* Symmetry **********************************

TEST_CASE("Symmetry mirrors velocity and passes scalars through", "[bc]")
{
    const auto symmetryVelocity = BoundaryTypeFactory::create
    (
        "symmetry", Field::Ux
    );
    const auto symmetryScalar = BoundaryTypeFactory::create
    (
        "symmetry", Field::p
    );

    REQUIRE(symmetryVelocity->typeName() == "symmetry");
    REQUIRE(symmetryVelocity->isSymmetry());
    REQUIRE(symmetryScalar->typeName() == "zeroGradient");
    REQUIRE(!symmetryScalar->isSymmetry());
}