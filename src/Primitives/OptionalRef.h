/******************************************************************************
 * @file OptionalRef.h
 * @brief Non-owning optional reference alias
 *
 * @details Defines OptionalRef<T> — a std::optional holding a
 * std::reference_wrapper<const T>. Used to express "may reference an
 * externally-owned T, or nothing" without raw pointers.
 *****************************************************************************/

#pragma once

#include <functional>
#include <optional>
#include <type_traits>


/// Non-owning optional reference (absent = std::nullopt)
template<typename T>
requires std::is_object_v<T>
using OptionalRef = std::optional<std::reference_wrapper<const T>>;
