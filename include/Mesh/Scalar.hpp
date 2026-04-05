/******************************************************************************
 * @file Scalar.hpp
 * @brief Floating-point precision configuration and tolerance definitions
 * 
 * @details This file defines the Scalar type used throughout the CFD solver 
 * and the numerical tolerances that adapt to the chosen precision.
 * Precision is controlled via CMake configuration.
 *****************************************************************************/

#pragma once

#include <string_view>
#include <limits>


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

/**
 * @brief Type-safe scalar literal conversion function
 * @tparam T Input type to convert
 * @param value Value to convert to Scalar type
 * @return Value converted to Scalar type
 * 
 * @details
 * This function ensures compile-time evaluation and type safety
 * when converting numeric literals to the configured Scalar type.
 */
template<typename T>
constexpr Scalar S(T value)
{    
    return static_cast<Scalar>(value);
}
