/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Field.cpp
 * @brief Implementation of Field identifier helpers
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Field.h"

// Project headers
#include "ErrorHandler.h"

// *************************** Non-Member Functions ***************************

NameRef fieldToString(Field field) noexcept
{
    using enum Field;
    switch (field)
    {
        case Ux:    return "Ux";
        case Uy:    return "Uy";
        case Uz:    return "Uz";
        case p:     return "p";
        case pCorr: return "pCorr";
        case k:     return "k";
        case omega: return "omega";
        case nut:   return "nut";
    }

    FatalError("Corrupted Field value");
}
