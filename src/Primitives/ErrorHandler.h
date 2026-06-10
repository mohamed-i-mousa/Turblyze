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
 *
 * @note FatalError is declared `[[noreturn]] noexcept` and calls
 * std::abort() unconditionally. It can therefore be invoked from any
 * function marked `noexcept` (including accessors and destructors) without
 * risking a silent std::terminate from a propagated exception.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <cstdlib>
#include <iostream>
#include <source_location>

#include "StringTypes.h"

// *********************************** Alias **********************************

using Location = std::source_location;

// ************************* Error Handling Functions *************************

/// Print a fatal error message and abort the program
[[noreturn]] inline void FatalError
(
    MessageRef errorMessage,
    const Location errorLocation = Location::current()
) noexcept
{
    std::cerr
        << '\n' << '\n' << "FATAL ERROR"
        << '\n' << "    " << errorLocation.file_name() << ':'
        << errorLocation.line()
        << '\n' << "    " << errorMessage << std::endl;

    std::abort();
}


/// Print a warning message and continue execution
inline void Warning
(
    MessageRef warningMessage,
    const Location warningLocation = Location::current()
) noexcept
{
    std::cerr
        << '\n' << "[WARNING] (" << warningLocation.file_name() << ':'
        << warningLocation.line() << ") " << warningMessage << std::endl;
}
