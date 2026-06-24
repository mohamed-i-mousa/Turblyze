/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file RuntimeSelection.h
 * @brief Helpers for validating and reporting runtime-selectable names
 *
 * @details Each scheme/solver/model family resolves a selected name in its own
 * create() and hand-writes the list of its selectable names (its available*()
 * accessor), which the case-file parser validates against. This header provides
 * the isKnown()/joinNames()/unknownSelection() helpers used to validate
 * selections and format name lists for diagnostics.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "StringTypes.h"

// ************************ namespace RuntimeSelection ************************

namespace RuntimeSelection
{

/// Whether selection matches one of the available names
[[nodiscard]] bool isKnown
(
    Name selection,
    const NameList& available
) noexcept;

/// Join names into a comma-separated list for messages
[[nodiscard]] Name joinNames(const NameList& available);

/// Abort with a "unknown <category> '<selection>'" message & listing available
[[noreturn]] void unknownSelection
(
    Name category,
    Name selection,
    const NameList& available
);

} // namespace RuntimeSelection
