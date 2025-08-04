#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <string>
#include <vector>
#include <map>

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
 * This function writes a VTK PolyData file using the VTK XML library.
 *
 * This function exports points and faces as polygons to a .vtp file
 * that can be visualized in ParaView. Cell-centered scalar fields are
 * mapped to face data for visualization.
 */
void writeVtkFile(const std::string& filename,
                  const std::vector<Vector>& allNodes,
                  const std::vector<Face>& allFaces,
                  const std::vector<Cell>& allCells,
                  const std::map<std::string, const ScalarField*>& scalarCellFields = {});

} // namespace VtkWriter

#endif // VTKWRITER_H
