/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryTypeFactory.h
 * @brief Runtime selection of boundary-condition types by case-file name
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <memory>

// Project headers
#include "BoundaryType.h"
#include "Field.h"
#include "Scalar.h"
#include "StringTypes.h"

// *************************** Forward Declarations ***************************

class CaseReader;

// *********************** namespace BoundaryTypeFactory **********************

namespace BoundaryTypeFactory
{

/// Create a boundary condition of the given type for the given field
[[nodiscard]] std::unique_ptr<BoundaryType> create
(
    const Name& typeName,
    Field field,
    Scalar value = S(0.0)
);

/// Case-file-selectable type names valid on the field
[[nodiscard]] NameList availableTypes(Field field);

} // namespace BoundaryTypeFactory
