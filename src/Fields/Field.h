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

#include <string_view>


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


/**
 * @brief Human-readable name of a field, for diagnostics and logging
 * @param field Field identifier
 * @return Field name
 * @note Terminates the program if the value is unrecognized
 * @note For diagnostics only, never use the returned string for dispatch
 */
[[nodiscard]] std::string_view fieldToString(Field field) noexcept;
