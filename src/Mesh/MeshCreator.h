/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshCreator.h
 * @brief Mesh read and geometry-preparation phase
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Project headers
#include "CaseConfiguration.h"
#include "Mesh.h"

// *************************** namespace MeshCreator **************************

namespace MeshCreator
{

/// Read, prepare, and optionally quality-check the configured mesh
[[nodiscard]] Mesh create(const CaseConfiguration& config);

} // namespace MeshCreator
