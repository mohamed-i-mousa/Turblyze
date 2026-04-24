/******************************************************************************
 * @file Scalar.h
 * @brief Floating-point precision configuration and tolerance definitions
 * 
 * @details This file defines the Scalar type used throughout the CFD solver 
 * and the numerical tolerances that adapt to the chosen precision.
 * Precision is controlled via CMake configuration.
 *****************************************************************************/

#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>


/// Floating-point precision type (configured via CMakeLists.txt)
#ifdef PROJECT_USE_DOUBLE_PRECISION
    using Scalar = double;
    constexpr std::string_view SCALAR_MODE = "double (FP64)";
#else
    using Scalar = float;
    constexpr std::string_view SCALAR_MODE = "float (FP32)";
#endif


/// Numerical tolerances
constexpr Scalar smallValue = std::numeric_limits<Scalar>::epsilon();

constexpr Scalar vSmallValue = std::numeric_limits<Scalar>::min();

constexpr Scalar largeValue = Scalar(1.0)/smallValue;


/// Numeric types that may be converted to Scalar via S()
template<typename T>
concept ScalarLiteral =
    std::floating_point<T>
 || std::same_as<T, int>
 || std::same_as<T, long>
 || std::same_as<T, long long>
 || std::same_as<T, int8_t>
 || std::same_as<T, size_t>;

/**
 * @brief Type-safe scalar literal conversion function
 * @tparam T Floating-point or integer literal type to convert
 * @param value Value to convert to Scalar type
 * @return Value converted to Scalar type
 */
template<ScalarLiteral T>
constexpr Scalar S(T value)
{
    return static_cast<Scalar>(value);
}
