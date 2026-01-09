/******************************************************************************
 * @file Scalar.hpp
 * @brief Floating-point precision configuration and tolerance definitions
 * 
 * This file defines the Scalar type used throughout the CFD solver,
 * along with numerical tolerances that adapt to the chosen precision.
 * Precision is controlled via CMake configuration.
 *****************************************************************************/

#ifndef SCALAR_HPP
#define SCALAR_HPP

#include <string>
#include <limits>

/// Floating-point precision type (configured via CMakeLists.txt)

#ifdef PROJECT_USE_DOUBLE_PRECISION
    using Scalar = double;
    constexpr std::string_view SCALAR_MODE = "double (FP64)";
#else
    using Scalar = float;
    constexpr std::string_view SCALAR_MODE = "float (FP32)";
#endif

// Numerical tolerances 

inline const Scalar smallValue =
    std::numeric_limits<Scalar>::epsilon();

inline const Scalar vSmallValue =
    std::numeric_limits<Scalar>::min();

// Minimum values for warnings 

inline const Scalar minArea = 1e-12;

inline const Scalar minVolume = 1e-30;

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
inline constexpr Scalar S(T value)
{    
    return static_cast<Scalar>(value);
}

#endif // SCALAR_HPP