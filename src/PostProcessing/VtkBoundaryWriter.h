/******************************************************************************
 * @file VtkBoundaryWriter.h
 * @brief VTK PolyData (.vtp) writer for wall boundary surfaces
 *
 * @details Exports wall boundary faces with face-centered scalar fields
 * (e.g. yPlus, wallShearStress) as a VTK PolyData surface suitable for
 * ParaView surface plots.
 *****************************************************************************/

#pragma once

#include <map>
#include <string>

#include "Scalar.h"
#include "Mesh.h"
#include "FaceData.h"


namespace VTK
{

/// Write wall boundary face data to VTK PolyData (.vtp) file
void writeWallBoundaryData
(
    const std::string& filename,
    const Mesh& mesh,
    const std::map<std::string,
    const FaceData<Scalar>*>& scalarFaceFields = {},
    bool debug = false
);

} // namespace VTK
