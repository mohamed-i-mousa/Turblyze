/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Field.h
 * @brief Solver field identifier
 *
 * @details Defines the Field enumeration used to identify solver fields
 * (velocity components, pressure, turbulence quantities) in boundary-condition
 * storage, BC lookups, and gradient reconstruction. It replaces error-prone
 * field-name string keys with a compiler-checked type.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "StringTypes.h"

// ***************************** enum class Field *****************************

enum class Field
{
    Ux,
    Uy,
    Uz,
    p,
    pCorr,
    k,
    omega,
    nut
};

// *************************** Non-Member Functions ***************************

/// Human-readable name of a field, for diagnostics and logging
[[nodiscard]] NameRef fieldToString(Field field) noexcept;
