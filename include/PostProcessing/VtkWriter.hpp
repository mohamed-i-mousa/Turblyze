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
 *
 * @see ParaView: https://www.paraview.org/
 * @see VTK Documentation: https://vtk.org/documentation/
 *****************************************************************************/

#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include <string>
#include <vector>
#include <map>
#include <set>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "CellData.hpp"

// VTK library integration requirement
typedef long long vtkIdType;


namespace VtkWriter
{

// Public API

/**
 * @brief Write simulation results to VTK UnstructuredGrid (.vtu) file
 *
 * @details This function exports 3D volumetric mesh and field data to
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
    const std::vector<Vector>& allNodes,
    const std::vector<Cell>& allCells,
    const std::vector<Face>& allFaces,
    const std::map<std::string,
    const ScalarField*>& scalarCellFields = {},
    const std::map<std::string,
    const VectorField*>& vectorCellFields = {},
    bool debug = false
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
ScalarField computeVorticityMagnitude
(
    const VectorField& vorticity
);

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
    const std::vector<Cell>& allCells
);

// Internal helper functions and types

/// Strain rate tensor components
struct StrainTensor
{
    Scalar S_11, S_22, S_33;  ///< Diagonal components
    Scalar S_12, S_13, S_23;  ///< Off-diagonal components

    /// Compute ||S||^2 = S_ij * S_ij
    Scalar normSquared() const
    {
        return S_11*S_11 + S_22*S_22 + S_33*S_33
             + 2.0 * (S_12*S_12 + S_13*S_13 + S_23*S_23);
    }
};

/**
 * @brief Generic helper to compute magnitude of any vector field
 * @param field Input vector field
 * @param fieldName Name for the output field
 * @return Scalar field containing magnitude at each cell
 */
ScalarField computeMagnitude
(
    const VectorField& field,
    const std::string& fieldName
);

/**
 * @brief Compute strain rate tensor from velocity gradient components
 * @param gradUx Gradient of x-velocity component at a cell
 * @param gradUy Gradient of y-velocity component at a cell
 * @param gradUz Gradient of z-velocity component at a cell
 * @return Strain rate tensor components
 */
StrainTensor computeStrainTensor
(
    const Vector& gradUx,
    const Vector& gradUy,
    const Vector& gradUz
);

/**
 * @brief Compute centroid of a face given node indices
 * @param faceNodes Node indices for the face
 * @param allNodes Vector of all node coordinates
 * @return Face centroid
 */
Vector computeFaceCentroid
(
    const std::vector<size_t>& faceNodes,
    const std::vector<Vector>& allNodes
);

/**
 * @brief Compute normal vector for a triangular face
 * @param triNodes Node indices for the triangle
 * @param allNodes Vector of all node coordinates
 * @return Face normal vector
 */
Vector computeTriangleNormal
(
    const std::vector<size_t>& triNodes,
    const std::vector<Vector>& allNodes
);

/**
 * @brief Compute normal vector for a quadrilateral face
 * @param quadNodes Node indices for the quadrilateral
 * @param allNodes Vector of all node coordinates
 * @return Face normal vector
 */
Vector computeQuadNormal
(
    const std::vector<size_t>& quadNodes,
    const std::vector<Vector>& allNodes
);

/**
 * @brief Order hexahedron nodes according to VTK convention
 *
 * @details VTK Hexahedron: nodes 0-3 form bottom quad, nodes 4-7 form
 * top quad. Uses face topology for robust ordering without
 * axis-alignment assumptions.
 *
 * @param faceNodeLists Node lists for all faces of the hexahedron
 * @param uniqueNodes Set of all unique node indices
 * @param allNodes Vector of all node coordinates
 * @return Ordered node indices for VTK, or empty vector if ordering fails
 */
std::vector<vtkIdType> orderHexahedronNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::set<size_t>& uniqueNodes,
    const std::vector<Vector>& allNodes
);

/**
 * @brief Order wedge (prism) nodes according to VTK convention
 *
 * @details VTK Wedge: nodes 0,1,2 form bottom triangle, nodes 3,4,5
 * form top triangle. Uses actual face topology to determine correct
 * node ordering.
 *
 * @param faceNodeLists Node lists for all faces of the wedge
 * @param uniqueNodes Set of all unique node indices
 * @param allNodes Vector of all node coordinates
 * @return Ordered node indices for VTK, or empty vector if ordering fails
 */
std::vector<vtkIdType> orderWedgeNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::set<size_t>& uniqueNodes,
    const std::vector<Vector>& allNodes
);

/**
 * @brief Order pyramid nodes according to VTK convention
 *
 * @details VTK Pyramid: nodes 0,1,2,3 form quad base, node 4 is apex.
 * Uses face topology for robust ordering without axis assumptions.
 *
 * @param faceNodeLists Node lists for all faces of the pyramid
 * @param uniqueNodes Set of all unique node indices
 * @param allNodes Vector of all node coordinates
 * @return Ordered node indices for VTK, or empty vector if ordering fails
 */
std::vector<vtkIdType> orderPyramidNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::set<size_t>& uniqueNodes,
    const std::vector<Vector>& allNodes
);

} // namespace VtkWriter

#endif // VTK_WRITER_HPP