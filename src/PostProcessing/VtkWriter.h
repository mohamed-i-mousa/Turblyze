#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <type_traits>

struct Vector;
struct Face;
struct Cell;

#include "Scalar.h" 
#include "Face.h"
#include "Cell.h"
#include "CellData.h"
#include "Vector.h"


namespace VtkWriter {

/**
 * This function writes a VTK legacy file for unstructured grids.
 *
 * This function exports points and cells to a .vtk file
 * that can be visualized in ParaView.
 *
 *  @param filename The path to the output .vtk file (e.g., "mesh_output.vtk").
 *  @param allNodes A constant reference to a vector of Point objects representing the mesh nodes.
 *  @param allFaces A constant reference to a vector of Face objects defining the connectivity of faces.
 *  @param allCells A constant reference to a vector of Cell objects defining the polyhedral cells.
 */
void writeVtkFile(const std::string& filename,
                  const std::vector<Vector>& allNodes,
                  const std::vector<Face>& allFaces,
                  const std::vector<Cell>& allCells,
                  const std::map<std::string, const ScalarField*>& scalarCellFields = {});

} // namespace VtkWriter

#endif // VTKWRITER_H
