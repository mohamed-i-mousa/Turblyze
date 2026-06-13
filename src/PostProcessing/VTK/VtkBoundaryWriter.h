/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file VtkBoundaryWriter.h
 * @brief VTK PolyData (.vtp) writer for boundary surfaces
 *
 * @details Exports all boundary faces with patch metadata and optional
 * face-centered scalar fields (e.g. yPlus, wallShearStress) as a VTK
 * PolyData surface suitable for ParaView surface plots.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library includes
#include <map>

// Project headers
#include "Scalar.h"
#include "Mesh.h"
#include "FaceData.h"
#include "StringTypes.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

/// Alias for face data maps
using FaceDataMap = std::map<Name, const FaceData<Scalar>*>;

/// Write boundary face data to VTK PolyData (.vtp) file
void writeBoundaryData
(
    const FilePath& filename,
    const Mesh& mesh,
    const FaceDataMap& scalarFaceFields = {},
    bool debug = false
);

} // namespace VTK
