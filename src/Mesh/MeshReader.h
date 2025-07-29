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
 * The function processes the following sections:
 * - Comments section (0): Skipped during parsing
 * - Dimension section (2): Determines if mesh is 2D or 3D (3D only supported)
 * - Nodes section (10): Reads node coordinates and populates allNodes vector
 * - Cells section (12): Reads cell declarations and allocates allCells vector
 * - Faces section (13): Reads face connectivity, owner/neighbor cells, and populates allFaces vector
 * - Boundaries section (45): Reads boundary patch information and populates allBoundaryPatches vector
 *
 * After reading the mesh file, the function:
 * 1. Populates cell data by establishing face-to-cell and cell-to-cell relationships
 * 2. Assigns face signs (+1 for owner cell, -1 for neighbor cell)
 * 3. Builds neighbor cell lists for each cell
 * 4. Validates mesh integrity (minimum faces per cell, minimum nodes per face)
 * 5. Prints summary statistics
 *
 * Parameters:
 * - filePath: Path to the Fluent mesh file (.msh format)
 * - allNodes: Output vector to store node coordinates (0-based indexing)
 * - allFaces: Output vector to store face connectivity and cell relationships
 * - allCells: Output vector to store cell data with face and neighbor indices
 * - allBoundaryPatches: Output vector to store boundary patch information
 *
 * The function clears all input vectors before populating them with new data.
 * All indices are converted to 0-based indexing internally.
 *
 * Throws std::runtime_error for:
 * - File opening failures
 * - Invalid hex/decimal conversions
 * - Index out of bounds errors
 * - Empty or malformed data lines
 * - 2D mesh files (not supported)
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
 */

void readMshFile(const std::string& filePath,
                 std::vector<Vector>& allNodes,
                 std::vector<Face>& allFaces,
                 std::vector<Cell>& allCells,
                 std::vector<BoundaryPatch>& allBoundaryPatches);

#endif 