/******************************************************************************
 * @file Field.h
 * @brief Solver field identifier
 *
 * @details Defines the Field enumeration used to identify solver fields
 * (velocity components, pressure, turbulence quantities) in boundary-condition
 * storage, BC lookups, and gradient reconstruction. It replaces error-prone
 * field-name string keys with a compiler-checked type.
 * 
 * @enum Field
 * - Enumeration of solver fields (Ux, Uy, Uz, p, pCorr, k, omega, nut)
 * ***************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <string_view>

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
[[nodiscard]] std::string_view fieldToString(Field field) noexcept;
