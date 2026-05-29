/******************************************************************************
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
#include <string>

// Project headers
#include "Scalar.h"
#include "Mesh.h"
#include "FaceData.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

/// Write boundary face data to VTK PolyData (.vtp) file
void writeBoundaryData
(
    const std::string& filename,
    const Mesh& mesh,
    const std::map<std::string,
    const FaceData<Scalar>*>& scalarFaceFields = {},
    bool debug = false
);

} // namespace VTK
