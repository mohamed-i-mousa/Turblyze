/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ErrorHandler.h
 * @brief Fatal error and warning functions for program diagnostics
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
