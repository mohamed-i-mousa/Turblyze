#ifndef SCALAR_H
#define SCALAR_H

#include <string>
#include <limits>

/**
 * @file Scalar.h
 * @brief Floating-point precision configuration and tolerance definitions
 * 
 * This file defines the Scalar type used throughout the CFD solver,
 * along with numerical tolerances that adapt to the chosen precision.
 * Precision is controlled via CMake configuration.
 */

/// Floating-point precision type (configured via CMakeLists.txt)
#ifdef PROJECT_USE_DOUBLE_PRECISION
    using Scalar = double;
    constexpr std::string_view SCALAR_MODE = "double (FP64)";
#else
    using Scalar = float;
    constexpr std::string_view SCALAR_MODE = "float (FP32)";
#endif

/// Tolerance for division by zero checks
inline const Scalar DIVISION_TOLERANCE = 
    std::numeric_limits<Scalar>::epsilon();

/// Tolerance for floating-point equality comparisons  
inline const Scalar EQUALITY_TOLERANCE = 
    std::numeric_limits<Scalar>::epsilon() * 100;

/// Tolerance for face area calculations
inline const Scalar AREA_TOLERANCE = 1e-12;

/// Tolerance for cell volume calculations
inline const Scalar VOLUME_TOLERANCE = 1e-30;

/// Tolerance for gradient calculations
inline const Scalar GRADIENT_TOLERANCE = 1e-12;

/**
 * @brief Type-safe scalar literal conversion function
 * @tparam T Input type to convert
 * @param val Value to convert to Scalar type
 * @return Value converted to Scalar type
 * 
 * This function ensures compile-time evaluation and type safety
 * when converting numeric literals to the configured Scalar type.
 */
template<typename T>
inline constexpr Scalar S(T val)
{    
    return static_cast<Scalar>(val);
}

#endif