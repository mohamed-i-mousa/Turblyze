/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryTypeFactory.cpp
 * @brief Runtime selection of boundary-condition types by case-file name
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "BoundaryTypeFactory.h"

// Project headers
#include "ErrorHandler.h"
#include "RuntimeSelection.h"
#include "FixedValue.h"
#include "FixedGradient.h"
#include "ZeroGradient.h"
#include "NoSlip.h"
#include "WallFunction.h"
#include "Symmetry.h"

// *********************** namespace BoundaryTypeFactory **********************

namespace BoundaryTypeFactory
{

std::unique_ptr<BoundaryType> create
(
    const Name& typeName,
    Field field,
    Scalar value
)
{
    if (!RuntimeSelection::isKnown(typeName, availableTypes(field)))
    {
        RuntimeSelection::unknownSelection
        (
            "boundary condition type for field '"
          + Name(fieldToString(field)) + "'",
            typeName,
            availableTypes(field)
        );
    }

    if (typeName == "fixedValue")
    {
        return std::make_unique<FixedValue>(value);
    }

    if (typeName == "noSlip")
    {
        return std::make_unique<NoSlip>();
    }

    if (typeName == "zeroGradient")
    {
        return std::make_unique<ZeroGradient>();
    }

    if (typeName == "fixedGradient")
    {
        return std::make_unique<FixedGradient>(value);
    }

    if (typeName == "symmetry")
    {
        if (field == Field::Ux || field == Field::Uy || field == Field::Uz)
        {
            return std::make_unique<Symmetry>(field);
        }
        return std::make_unique<ZeroGradient>();
    }

    // The remaining selectable tokens are the wall-function flavors
    return std::make_unique<WallFunction>(typeName);
}


NameList availableTypes(Field field)
{
    switch (field)
    {
        case Field::Ux:
        case Field::Uy:
        case Field::Uz:
            return
                {
                    "fixedGradient",
                    "fixedValue",
                    "noSlip",
                    "symmetry",
                    "zeroGradient"
                };

        case Field::p:
            return
                {
                    "fixedGradient",
                    "fixedValue",
                    "symmetry",
                    "zeroGradient"
                };

        case Field::k:
            return
                {
                    "fixedGradient",
                    "fixedValue",
                    "kWallFunction",
                    "symmetry",
                    "zeroGradient"
                };

        case Field::omega:
            return
                {
                    "fixedGradient",
                    "fixedValue",
                    "omegaWallFunction",
                    "symmetry",
                    "zeroGradient"
                };

        case Field::nut:
            return
                {
                    "fixedGradient",
                    "fixedValue",
                    "nutWallFunction",
                    "symmetry",
                    "zeroGradient"
                };

        case Field::pCorr:
        default:
            return {};
    }
}

} // namespace BoundaryTypeFactory
