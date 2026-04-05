/******************************************************************************
 * @file ErrorHandler.hpp
 * @brief Fatal error and warning macros for program diagnostics
 *
 * Provides two macros for error reporting:
 *
 * - FatalError("message")  — prints file, line, and message to stderr,
 *                            then calls std::abort() (enables core dumps)
 * - Warning("message")     — prints file, line, and message to stderr,
 *                            then continues execution
 *
 * CFD solvers are batch executables, not servers. There is no caller to
 * catch and retry, so unrecoverable errors terminate immediately rather
 * than throwing exceptions. std::abort() preserves the full call stack
 * for post-mortem debugging.
 *****************************************************************************/

#pragma once

#include <cstdlib>
#include <iostream>
#include <string>


/// Print a fatal error message and abort the program
[[noreturn]] inline void fatalErrorImpl
(
    const char* file,
    int line,
    const std::string& msg
)
{
    std::cerr
        << "\n\nFATAL ERROR"
        << "\n    " << file << ':' << line
        << "\n    " << msg
        << '\n' << std::endl;

    std::abort();
}


/// Print a warning message and continue execution
inline void warningImpl
(
    const char* file,
    int line,
    const std::string& msg
)
{
    std::cerr
        << "\n[WARNING] (" << file << ':' << line << ") "
        << msg << std::endl;
}


// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define FatalError(msg) fatalErrorImpl(__FILE__, __LINE__, msg)

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define Warning(msg) warningImpl(__FILE__, __LINE__, msg)
