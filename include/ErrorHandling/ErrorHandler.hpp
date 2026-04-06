/******************************************************************************
 * @file ErrorHandler.hpp
 * @brief Fatal error and warning functions for program diagnostics
 *
 * Provides two functions for error reporting:
 *
 * - FatalError("message")  — prints file, line, and message to stderr,
 *                            then calls std::abort() (enables core dumps)
 * - Warning("message")     — prints file, line, and message to stderr,
 *                            then continues execution
 *
 * The CFD solver is batch executable, not a server. There is no caller to
 * catch and retry, so unrecoverable errors terminate immediately rather
 * than throwing exceptions. std::abort() preserves the entire memory snapshot
 * for post-crash debugging.
 *****************************************************************************/

#pragma once

#include <cstdlib>
#include <iostream>
#include <source_location>
#include <string>


/**
 * @brief Print a fatal error message and abort the program
 * @param message The error message to display
 * @param location The source location (file & line) where the error occurred
 */
// Print a fatal error message and abort the program
[[noreturn]] inline void FatalError
(
    const std::string& message,
    std::source_location location = std::source_location::current()
)
{
    std::cerr
        << "\n\nFATAL ERROR"
        << "\n    " << location.file_name() << ':' << location.line()
        << "\n    " << message
        << '\n' << std::endl;

    std::abort();
}


/**
 * @brief Print a warning message and continue execution
 * @param message The warning message to display
 * @param location The source location (file & line) where the warning occurred
 */
inline void Warning
(
    const std::string& message,
    std::source_location location = std::source_location::current()
)
{
    std::cerr
        << "\n[WARNING] (" << location.file_name() << ':' 
        << location.line() << ") " << message << std::endl;
}
