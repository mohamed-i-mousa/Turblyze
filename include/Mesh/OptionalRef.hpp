/******************************************************************************
 * @file OptionalRef.hpp
 * @brief Non-owning optional reference alias
 *
 * @details Defines OptionalRef<T> — a std::optional holding a
 * std::reference_wrapper<const T>. Used to express "may reference an
 * externally-owned T, or nothing" without raw pointers.
 *****************************************************************************/

#pragma once

#include <optional>
#include <functional>


/// Non-owning optional reference (absent = std::nullopt)
template<typename T>
using OptionalRef = std::optional<std::reference_wrapper<const T>>;
