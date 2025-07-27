#ifndef MESHREADER_H
#define MESHREADER_H

#include <string>
#include <vector>
#include "Scalar.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"

/*
 * This function reads a Fluent mesh file exported from ANSYS Meshing.
 *
 * The mesh file is expected to follow a specific format with sections for nodes,
 * faces, and cells. Numerical values like IDs, counts, and connectivity are
 * expected in hexadecimal, while point coordinates are in decimal.
 *
 *  In faces section:       
 *       Type:
 *       "2" = "internal"
 *       "3" = "wall";
 *       "4" = "Pressure-Inlet / Inlet-Vet / Intake-Fan"
 *       "5" = "Pressure-Outet / Exhaust-Fan / Outlet-Vent"
 *       "7" = "Symmetry"
 *       "8" = "Periodic-Shadow"
 *       "9" = "Pressure-Far-Field"
 *       "a" = "Velocity-Inlet"
 *       "c" = "Periodic"
 *       "e" = "Fan / Porous-Jump / Radiator"
 *      "14" = "Mass-Flow-Inlet"
 *      "18" = "Interface"
 *      "1F" = "Parent"
 *      "24" = "Outflow"
 *      "25" = "Axis"
 *
 *       Element Type:
 *       "0" = "mixed"
 *       "2" = "Line/Edge"
 *       "3" = "Triangular"
 *       "4" = "Quadrilateral"
 *       "5" = "Polygonal"
 *       
 *
 * filePath Path to the mesh file.
 * allNodes Vector to be filled with Scalar data.
 * allFaces Vector to be filled with Face data.
 * allCells Vector to be filled with Cell data.
 */

void readMshFile(const std::string& filePath,
                 std::vector<Vector>& allNodes,
                 std::vector<Face>& allFaces,
                 std::vector<Cell>& allCells,
                 std::vector<BoundaryPatch>& allBoundaryPatches);

#endif 