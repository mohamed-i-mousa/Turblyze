/******************************************************************************
 * @file ErrorHandler.h
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
#include <string_view>


/// Print a fatal error message and abort the program
[[noreturn]] inline void FatalError
(
    std::string_view message,
    std::source_location location = std::source_location::current()
) noexcept
{
    std::cerr
        << '\n' << '\n' << "FATAL ERROR"
        << '\n' << "    " << location.file_name() << ':' << location.line()
        << '\n' << "    " << message
        << std::endl;

    std::abort();
}


/// Print a warning message and continue execution
inline void Warning
(
    std::string_view message,
    std::source_location location = std::source_location::current()
) noexcept
{
    std::cerr
        << '\n' << "[WARNING] (" << location.file_name() << ':'
        << location.line() << ") " << message << std::endl;
}
