/******************************************************************************
 * @file VtkWriter.h
 * @brief VTK UnstructuredGrid (.vtu) writer for 3D volume cell data
 *
 * @details Exports the 3D unstructured mesh (tetrahedra, hexahedra, wedges,
 * pyramids) together with cell-centered scalar and vector fields in the
 * VTK UnstructuredGrid format, suitable for ParaView volume rendering,
 * slicing, clipping, and isosurface extraction.
 *
 * @see ParaView: https://www.paraview.org/
 * @see VTK Documentation: https://vtk.org/documentation/
 *****************************************************************************/

#pragma once

#include <array>
#include <map>
#include <string>

#include "Scalar.h"
#include "Mesh.h"
#include "CellData.h"


namespace VTK
{

/// Write simulation results to VTK UnstructuredGrid (.vtu) file
void writeVtkUnstructuredGrid
(
    const std::string& filename,
    const Mesh& mesh,
    const std::map<std::string,
    const ScalarField*>& scalarCellFields = {},
    const std::map<std::string,
    std::array<const ScalarField*, 3>>& vectorCellFields = {},
    bool debug = false
);

} // namespace VTK
