/******************************************************************************
 * @file VtkWriter.hpp
 * @brief VTK (Visualization Toolkit) File Writer
 *
 * @details This header defines the VtkWriter namespace, a collection of
 * utility functions for exporting simulation results in VTK format compatible
 * with ParaView.
 *
 * The VtkWriter namespace provides:
 * - VTK UnstructuredGrid export (.vtu) for 3D volume cells (tetrahedra,
 *   hexahedra, wedges, pyramids)
 * - Cell-centered scalar fields (pressure, turbulence quantities)
 * - Cell-centered vector fields (velocity, gradients)
 * - PVD time series files for transient animations
 * - Derived field computation (velocity/vorticity magnitude, Q-criterion,
 *   strain rate magnitude)
 *
 * @see ParaView: https://www.paraview.org/
 * @see VTK Documentation: https://vtk.org/documentation/
 *****************************************************************************/

#pragma once

#include <string>
#include <vector>
#include <map>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"


namespace VtkWriter
{

/**
 * @brief Write simulation results to VTK UnstructuredGrid (.vtu) file
 *
 * @details
 * This function exports 3D volumetric mesh and field data to
 * VTK UnstructuredGrid format (.vtu) which enables full 3D visualization
 * including volume rendering, slicing, clipping, and isosurfaces.
 *
 * UnstructuredGrid exports the actual 3D cells (tetrahedra, hexahedra,
 * prisms, pyramids) with proper cell-centered data.
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
    std::span<const Vector> allNodes,
    std::span<const Cell> allCells,
    std::span<const Face> allFaces,
    const std::map<std::string,
    const ScalarField*>& scalarCellFields = {},
    const std::map<std::string,
    const VectorField*>& vectorCellFields = {},
    bool debug = false
);

/**
 * @brief Write wall boundary face data to VTK PolyData (.vtp) file
 *
 * @details
 * Exports wall boundary faces with face-centered scalar fields (e.g.
 * yPlus, wallShearStress) as VTK PolyData for surface visualization.
 *
 * @param filename Output VTK file path (should end with .vtp extension)
 * @param allNodes Vector of 3D node coordinates
 * @param allFaces Vector of mesh faces (wall faces are extracted)
 * @param scalarFaceFields Map of scalar field names to face-centered data
 * @param debug Enable verbose output
 */
void writeWallBoundaryData
(
    const std::string& filename,
    std::span<const Vector> allNodes,
    std::span<const Face> allFaces,
    const std::map<std::string,
    const FaceData<Scalar>*>& scalarFaceFields = {},
    bool debug = false
);

/**
 * @brief Compute velocity magnitude field from velocity vector field
 * @param velocity Input velocity vector field
 * @return Scalar field containing velocity magnitude at each cell
 */
[[nodiscard("Computed derived field needed for visualization output")]]
ScalarField computeVelocityMagnitude(const VectorField& velocity);

/**
 * @brief Compute vorticity magnitude field from vorticity vector field
 * @param vorticity Input vorticity vector field
 * @return Scalar field containing vorticity magnitude at each cell
 */
[[nodiscard("Computed derived field needed for visualization output")]]
ScalarField computeVorticityMagnitude
(
    const VectorField& vorticity
);

/**
 * @brief Compute Q-criterion for vortex identification
 *
 * @details
 * Q-criterion identifies vortex cores as regions where Q > 0,
 * where Q = 0.5 * (||Omega||^2 - ||S||^2), with Omega being the
 * rotation rate tensor and S the strain rate tensor.
 *
 * @param gradUx Gradient of x-velocity component
 * @param gradUy Gradient of y-velocity component
 * @param gradUz Gradient of z-velocity component
 * @return Scalar field containing Q-criterion at each cell
 */
[[nodiscard("Computed derived field needed for visualization output")]]
ScalarField computeQCriterion
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
);

/**
 * @brief Compute strain rate magnitude field
 *
 * @details
 * Strain rate magnitude = sqrt(2 * S_ij * S_ij) where
 * S_ij is the symmetric strain rate tensor.
 *
 * @param gradUx Gradient of x-velocity component
 * @param gradUy Gradient of y-velocity component
 * @param gradUz Gradient of z-velocity component
 * @return Scalar field containing strain rate magnitude at each cell
 */
[[nodiscard("Computed derived field needed for visualization output")]]
ScalarField computeStrainRateMagnitude
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
);

/**
 * @brief Create PVD time series file header for transient runs
 * @param pvdFilename Path to .pvd file to create
 */
void writePVDTimeSeriesHeader(const std::string& pvdFilename);

/**
 * @brief Append a timestep to PVD time series file
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

/**
 * @brief Write cell geometry data (volumes and centroids) to text file
 * @param filename Output text file path
 * @param allCells Vector of mesh cells with computed geometry
 */
void writeCellGeometryData
(
    const std::string& filename,
    std::span<const Cell> allCells
);

} // namespace VtkWriter
