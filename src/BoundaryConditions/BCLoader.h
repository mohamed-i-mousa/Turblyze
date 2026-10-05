/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BCLoader.h
 * @brief Case-file boundary condition registration
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "BoundaryConditions.h"
#include "CaseConfiguration.h"
#include "Mesh.h"

// *************************** Forward Declarations ***************************

class CaseReader;

// **************************** namespace BCLoader ****************************

namespace BCLoader
{

/// Parse and construct the boundary conditions from the case file
[[nodiscard]] BoundaryConditions load
(
    const CaseReader& reader,
    const CaseConfiguration& config,
    const Mesh& mesh
);

} // namespace BCLoader
