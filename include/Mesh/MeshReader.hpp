
/******************************************************************************
 * @file MeshReader.h
 * @brief Fluent mesh file reader for ANSYS mesh files
 *****************************************************************************/

#ifndef MESHREADER_H
#define MESHREADER_H

#include <string>
#include <vector>
#include "Scalar.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "BoundaryPatch.hpp"

/**
 * @brief Read Fluent mesh file exported from ANSYS Meshing
 * @param filePath Path to the Fluent mesh file (.msh format)
 * @param allNodes Output vector to store node coordinates (0-based indexing)
 * @param allFaces Output vector to store face connectivity and relationships
 * @param allCells Output vector to store cell data with connectivity
 * @param allBoundaryPatches Output vector to store boundary patch information
 * 
 * @throws std::runtime_error File opening failures
 * @throws std::runtime_error Invalid hex/decimal conversions
 * @throws std::runtime_error Index out of bounds errors
 * @throws std::runtime_error Empty or malformed data lines
 * @throws std::runtime_error 2D mesh files (not supported)
 * 
 * This function reads Fluent mesh files with the following sections:
 * - Comments section (0): Skipped during parsing
 * - Dimension section (2): Determines if mesh is 2D or 3D (3D only)
 * - Nodes section (10): Reads node coordinates
 * - Cells section (12): Reads cell declarations and allocates storage
 * - Faces section (13): Reads face connectivity and owner/neighbor cells
 * - Boundaries section (45): Reads boundary patch information
 * 
 * Post-processing steps:
 * 1. Establishes face-to-cell and cell-to-cell relationships
 * 2. Assigns face signs (+1 for owner cell, -1 for neighbor cell)
 * 3. Builds neighbor cell lists for each cell
 * 4. Validates mesh integrity
 * 5. Prints summary statistics
 * 
 * All indices are converted to 0-based indexing internally.
 * Numerical values (IDs, counts, connectivity) are in hexadecimal format.
 * Point coordinates are in decimal format.
 * 
 * @note Fluent face types (hexadecimal):
 * - "2" = internal, "3" = wall, "4" = pressure-inlet
 * - "5" = pressure-outlet, "7" = symmetry, "8" = periodic-shadow
 * - "9" = pressure-far-field, "a" = velocity-inlet, "c" = periodic
 * - "e" = fan/porous-jump, "14" = mass-flow-inlet, "18" = interface
 * - "1F" = parent, "24" = outflow, "25" = axis
 * 
 * @note Fluent element types (hexadecimal):
 * - "0" = mixed, "2" = line/edge, "3" = triangular
 * - "4" = quadrilateral, "5" = polygonal
 */
void readMshFile
(
    const std::string& filePath,
    std::vector<Vector>& allNodes,
    std::vector<Face>& allFaces,
    std::vector<Cell>& allCells,
    std::vector<BoundaryPatch>& allBoundaryPatches
);

#endif