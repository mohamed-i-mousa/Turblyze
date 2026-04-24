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

#include <map>
#include <string>

#include "Scalar.h"
#include "Mesh.h"
#include "CellData.h"


namespace VTK
{

/**
 * @brief Write simulation results to VTK UnstructuredGrid (.vtu) file
 *
 * @details Exports 3D volumetric mesh and field data to VTK
 * UnstructuredGrid format (.vtu) which enables full 3D visualization
 * including volume rendering, slicing, clipping, and isosurfaces.
 *
 * @param filename Output VTK file path (should end with .vtu extension)
 * @param mesh Mesh view (nodes, faces, cells)
 * @param scalarCellFields Map of scalar field names to cell-centered data
 * @param vectorCellFields Map of vector field names to cell-centered data
 * @param debug Enable verbose output of cell-type diagnostics
 */
void writeVtkUnstructuredGrid
(
    const std::string& filename,
    const Mesh& mesh,
    const std::map<std::string,
    const ScalarField*>& scalarCellFields = {},
    const std::map<std::string,
    const VectorField*>& vectorCellFields = {},
    bool debug = false
);

} // namespace VTK
