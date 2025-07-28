#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>

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
 */
void writeVtkFile(const std::string& filename,
                  const std::vector<Vector>& allNodes,
                  const std::vector<Face>& allFaces,
                  const std::vector<Cell>& allCells,
                  const std::map<std::string, const ScalarField*>& scalarCellFields = {});

} // namespace VtkWriter

#endif // VTKWRITER_H
