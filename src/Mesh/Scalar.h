#ifndef SCALAR_H
#define SCALAR_H

#include <string>
#include <limits>

// Floating-point precision -- Defined in CMakeLists.txt
#ifdef PROJECT_USE_DOUBLE_PRECISION
    using Scalar = double;
    constexpr std::string_view SCALAR_MODE = "double (FP64)";
#else
    using Scalar = float;
    constexpr std::string_view SCALAR_MODE = "float (FP32)";
#endif

// Global numerical tolerances - adapt automatically to Scalar precision
inline const Scalar DIVISION_TOLERANCE = std::numeric_limits<Scalar>::epsilon();
inline const Scalar EQUALITY_TOLERANCE = std::numeric_limits<Scalar>::epsilon() * 100;
inline const Scalar AREA_TOLERANCE = 1e-12;
inline const Scalar VOLUME_TOLERANCE = 1e-30;
inline const Scalar GRADIENT_TOLERANCE = 1e-12;

// Literal conversion
template<typename T>                  // template is used to allow the function to be used with different types
inline constexpr Scalar S(T val) {    // constexpr is used to ensure that the function is evaluated at compile time
    return static_cast<Scalar>(val);
}

#endif