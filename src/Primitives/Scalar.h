/******************************************************************************
 * @file Scalar.h
 * @brief Floating-point precision configuration and tolerance definitions
 *
 * @details This file defines the Scalar type used throughout the CFD solver
 * and the numerical tolerances that adapt to the chosen precision.
 * Precision is controlled via CMake configuration.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <concepts>
#include <cstdint>
#include <limits>
#include <vector>

// Project headers
#include "Integer.h"
#include "StringTypes.h"

// ************************** Precision Configuration *************************

/// Floating-point precision type (configured via CMakeLists.txt)
#ifdef PROJECT_USE_DOUBLE_PRECISION
    using Scalar = double;
    constexpr MessageRef SCALAR_MODE = "double (FP64)";
#else
    using Scalar = float;
    constexpr MessageRef SCALAR_MODE = "float (FP32)";
#endif

// ***************************** Scalar Conversion ****************************

/// Numeric types that may be converted to Scalar via S()
template<typename T>
concept ScalarLiteral =
    std::floating_point<T>
 || std::same_as<T, int>
 || std::same_as<T, long>
 || std::same_as<T, long long>
 || std::same_as<T, std::int8_t>
 || std::same_as<T, Count>;

/// Type-safe scalar literal conversion function
template<ScalarLiteral T>
[[nodiscard]] constexpr Scalar S(T value) noexcept
{
    return static_cast<Scalar>(value);
}

// *************************** Numerical Tolerances ***************************

constexpr Scalar smallValue = std::numeric_limits<Scalar>::epsilon();

constexpr Scalar vSmallValue = std::numeric_limits<Scalar>::min();

constexpr Scalar largeValue = S(1.0) / smallValue;

// ********************************* Aliases **********************************

using ScalarList = std::vector<Scalar>;
