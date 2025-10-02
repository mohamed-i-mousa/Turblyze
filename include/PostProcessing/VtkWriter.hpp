/******************************************************************************
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
 *****************************************************************************/

#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <string>
#include <vector>
#include <map>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "CellData.hpp"

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
    const ScalarField*>& scalarCellFields = {},
    const std::map<std::string,
    const VectorField*>& vectorCellFields = {}
);

/**
 * @brief Write CFD simulation results to VTK UnstructuredGrid file
 *
 * @details This function exports 3D volumetric mesh and field data to VTK
 * UnstructuredGrid format (.vtu) which enables full 3D visualization including
 * volume rendering, slicing, clipping, and isosurfaces.
 *
 * Unlike PolyData which only represents surfaces, UnstructuredGrid exports
 * the actual 3D cells (tetrahedra, hexahedra, prisms, pyramids) with proper
 * cell-centered data. This is the recommended format for CFD volume visualization.
 *
 * @param filename Output VTK file path (should end with .vtu extension)
 * @param allNodes Vector of 3D node coordinates
 * @param allCells Vector of mesh cells with connectivity
 * @param allFaces Vector of mesh faces (needed for cell node extraction)
 * @param scalarCellFields Map of scalar field names to cell-centered data
 * @param vectorCellFields Map of vector field names to cell-centered data
 */
void writeVtkUnstructuredGrid
(
    const std::string& filename,
    const std::vector<Vector>& allNodes,
    const std::vector<Cell>& allCells,
    const std::vector<Face>& allFaces,
    const std::map<std::string,
    const ScalarField*>& scalarCellFields = {},
    const std::map<std::string,
    const VectorField*>& vectorCellFields = {}
);

/**
 * @brief Compute velocity magnitude field from velocity vector field
 * @param velocity Input velocity vector field
 * @return Scalar field containing velocity magnitude at each cell
 */
ScalarField computeVelocityMagnitude(const VectorField& velocity);

/**
 * @brief Compute vorticity magnitude field from vorticity vector field
 * @param vorticity Input vorticity vector field
 * @return Scalar field containing vorticity magnitude at each cell
 */
ScalarField computeVorticityMagnitude(const VectorField& vorticity);

/**
 * @brief Compute Q-criterion for vortex identification
 *
 * @details Q-criterion identifies vortex cores as regions where Q > 0,
 * where Q = 0.5 * (||Omega||^2 - ||S||^2), with Omega being the
 * rotation rate tensor and S the strain rate tensor.
 *
 * @param gradUx Gradient of x-velocity component
 * @param gradUy Gradient of y-velocity component
 * @param gradUz Gradient of z-velocity component
 * @return Scalar field containing Q-criterion at each cell
 */
ScalarField computeQCriterion
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
);

/**
 * @brief Compute strain rate magnitude field
 *
 * @details Strain rate magnitude = sqrt(2 * S_ij * S_ij) where
 * S_ij is the symmetric strain rate tensor.
 *
 * @param gradUx Gradient of x-velocity component
 * @param gradUy Gradient of y-velocity component
 * @param gradUz Gradient of z-velocity component
 * @return Scalar field containing strain rate magnitude at each cell
 */
ScalarField computeStrainRateMagnitude
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
);

/**
 * @brief Create PVD time series file header
 *
 * @details Creates a ParaView Data (PVD) file for time series animations.
 * Call this once at the beginning, then use appendPVDTimeStep() to add
 * timesteps.
 *
 * @param pvdFilename Path to .pvd file to create
 */
void writePVDTimeSeriesHeader(const std::string& pvdFilename);

/**
 * @brief Append a timestep to PVD time series file
 *
 * @details Adds a new timestep entry to an existing PVD file. The VTU
 * file should already be written before calling this function.
 *
 * @param pvdFilename Path to existing .pvd file
 * @param vtuFilename Relative path to .vtu file for this timestep
 * @param timeValue Physical time value for this timestep
 */
void appendPVDTimeStep
(
    const std::string& pvdFilename,
    const std::string& vtuFilename,
    Scalar timeValue
);

} // namespace VtkWriter

#endif // VTKWRITER_H
