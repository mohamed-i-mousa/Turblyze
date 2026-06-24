/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file RuntimeSelection.cpp
 * @brief Queries over the runtime-selectable component name lists
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "RuntimeSelection.h"

// Project headers
#include "ErrorHandler.h"

// *********************** namespace RuntimeSelection ************************

namespace RuntimeSelection
{

bool isKnown
(
    Name selection,
    const NameList& available
) noexcept
{
    for (const Name& name : available)
    {
        if (selection == name)
        {
            return true;
        }
    }

    return false;
}


Name joinNames(const NameList& available)
{
    Name list;
    for (const Name& name : available)
    {
        if (!list.empty())
        {
            list += ", ";
        }
        list += name;
    }

    return list;
}


void unknownSelection
(
    Name category,
    Name selection,
    const NameList& available
)
{
    FatalError
    (
        "Unknown " + category + " '" + selection
      + "'. Valid options: " + joinNames(available) + "."
    );
}

} // namespace RuntimeSelection
