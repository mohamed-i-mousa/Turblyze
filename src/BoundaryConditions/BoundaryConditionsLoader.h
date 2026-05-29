/******************************************************************************
 * @file BoundaryConditionLoader.h
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

/// Register all mesh patches and configured field boundary conditions
void load
(
    const CaseReader& reader,
    const CaseConfiguration& config,
    Mesh& mesh,
    BoundaryConditions& bcManager
);

} // namespace BCLoader
