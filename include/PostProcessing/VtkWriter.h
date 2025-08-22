/**
 * @file VtkWriter.h
 * @brief VTK (Visualization Toolkit) File Writer for CFD Results
 * 
 * @details This header file defines the VtkWriter namespace, which provides 
 * functionality to export CFD simulation results in VTK format for 
 * visualization in ParaView and other VTK-compatible visualization tools.
 * 
 * The VtkWriter supports:
 * - Export of mesh geometry (nodes and faces)
 * - Cell-centered scalar fields (pressure, temperature, turbulence quantities)
 * - Vector fields (velocity, gradients)
 * - Proper data mapping from cell centers to face centers for visualization
 * - VTK PolyData format (.vtp files) for efficient rendering
 * 
 * @see ParaView: https://www.paraview.org/
 * @see VTK Documentation: https://vtk.org/documentation/
 * @see VTK File Formats: https://vtk.org/Wiki/VTK_XML_Formats
 */

#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <string>
#include <vector>
#include <map>

#include "Scalar.h" 
#include "Vector.h"
#include "Face.h"
#include "CellData.h"

namespace VtkWriter {

/**
 * @brief Write CFD simulation results to VTK PolyData file
 * 
 * @details This function exports the complete CFD simulation results including
 * mesh geometry and field data to a VTK PolyData file (.vtp) that can be
 * directly opened in ParaView for visualization and analysis.
 * 
 * The function performs the following operations:
 * 1. **Mesh Export**: Writes node coordinates and face connectivity
 * 2. **Data Mapping**: Maps cell-centered fields to face data for visualization
 * 3. **Field Export**: Exports scalar fields (pressure, temperature, etc.)
 * 4. **File Generation**: Creates a valid VTK XML file with proper formatting
 * 
 * **Data Mapping Strategy**:
 * Since VTK PolyData represents surfaces using faces, cell-centered data
 * must be mapped to face centers. The mapping uses:
 * - Linear interpolation between adjacent cell centers
 * - Proper handling of boundary faces
 * - Conservation of field values during mapping
 * 
 * **Supported Field Types**:
 * - Scalar fields: pressure, temperature, turbulent kinetic energy
 * - Vector fields: velocity, gradients, forces
 * - Derived quantities: vorticity, strain rate, etc.
 * 
 * **File Format**:
 * The output file follows VTK XML PolyData format:
 * ```xml
 * <?xml version="1.0"?>
 * <VTKFile type="PolyData" version="0.1">
 *   <PolyData>
 *     <Piece NumberOfPoints="N" NumberOfPolys="M">
 *       <Points>...</Points>
 *       <Polys>...</Polys>
 *       <PointData>...</PointData>
 *       <CellData>...</CellData>
 *     </Piece>
 *   </PolyData>
 * </VTKFile>
 * ```
 * 
 * @param filename Output VTK file path (should end with .vtp extension)
 * @param allNodes Vector of 3D node coordinates defining the mesh
 * @param allFaces Vector of mesh faces with connectivity information
 * @param scalarCellFields Map of field names to cell-centered scalar fields
 * 
 * @note The filename should have a .vtp extension for proper ParaView 
 * recognition
 * @note Cell-centered fields are automatically mapped to face centers for 
 * visualization
 * @note Empty scalarCellFields map results in geometry-only export
 * 
 * @example
 * ```cpp
 * // Export mesh with pressure and velocity fields
 * std::map<std::string, const ScalarField*> fields;
 * fields["pressure"] = &pressureField;
 * fields["velocity_magnitude"] = &velocityMagnitudeField;
 * 
 * VtkWriter::writeVtkFile("results.vtp", nodes, faces, fields);
 * ```
 * 
 * @throws std::runtime_error if file cannot be created or written
 * @throws std::invalid_argument if mesh data is invalid
 * 
 * @return void
 * 
 * @see writeVtkFileUnstructured() for unstructured grid export
 * @see writeVtkFileStructured() for structured grid export
 * @see exportFieldData() for field-only export
 */
void writeVtkFile
(
    const std::string& filename,
    const std::vector<Vector>& allNodes,
    const std::vector<Face>& allFaces,
    const std::map<std::string, 
    const ScalarField*>& scalarCellFields = {}
);

} // namespace VtkWriter

#endif // VTKWRITER_H
