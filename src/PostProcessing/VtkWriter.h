#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <string>
#include <vector>
#include <map>

#include "Scalar.h" 
#include "Vector.h"
#include "Face.h"
#include "Vector.h"
#include "CellData.h"

namespace VtkWriter {

/**
 * This function writes a VTK PolyData file using the VTK XML library.
 *
 * This function exports points and faces as polygons to a .vtp file
 * that can be visualized in ParaView. Cell-centered scalar fields are
 * mapped to face data for visualization.
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
