/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshReader.h
 * @brief Fluent mesh file reader for ANSYS mesh files
 *
 * @details MeshReader reads mesh data from Fluent (.msh) files and converts
 * them into the internal structure (Nodes, Faces, Cells). Currently supports
 * 3D unstructured meshes exported from ANSYS Meshing.
 *
 * @note Supported Fluent face types (hexadecimal):
 * - "2" = internal, "3" = wall, "4" = pressure-inlet
 * - "5" = pressure-outlet, "7" = symmetry, "8" = periodic-shadow
 * - "9" = pressure-far-field, "a" = velocity-inlet, "c" = periodic
 * - "e" = fan/porous-jump, "14" = mass-flow-inlet, "18" = interface
 * - "1F" = parent, "24" = outflow, "25" = axis
 *
 * @note Supported Fluent element types (hexadecimal):
 * - "0" = mixed, "2" = line/edge, "3" = triangular
 * - "4" = quadrilateral, "5" = polygonal
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <array>
#include <iosfwd>
#include <utility>

// Project headers
#include "MeshContainers.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"
#include "Integer.h"
#include "StringTypes.h"

// ***************************** class MeshReader *****************************

class MeshReader
{
public:

    using TokenList = std::vector<Token>;

// ************************* Special Member Functions *************************

    /// Construct MeshReader and parse the given Fluent mesh file
    explicit MeshReader(const FilePath& filePath);

// ***************************** Accessor Methods *****************************

    /// Transfer ownership of nodes data
    [[nodiscard]] NodeList moveNodes() noexcept
    {
        return std::move(nodes_);
    }

    /// Transfer ownership of faces data
    [[nodiscard]] FaceList moveFaces() noexcept
    {
        return std::move(faces_);
    }

    /// Transfer ownership of cells data
    [[nodiscard]] CellList moveCells() noexcept
    {
        return std::move(cells_);
    }

    /// Transfer ownership of boundary patches data
    [[nodiscard]] PatchList moveBoundaryPatches() noexcept
    {
        return std::move(boundaryPatches_);
    }

// ****************************** Private Members *****************************

private:

    /// All mesh node coordinates
    NodeList nodes_;

    /// All mesh faces
    FaceList faces_;

    /// All mesh cells
    CellList cells_;

    /// All boundary patches
    PatchList boundaryPatches_;

// ****************************** Private Methods *****************************

private:

    /// Parse the complete mesh file
    void parseFile(const FilePath& filePath);

    /// Parse the comment section and skip its contents
    void parseCommentSection(std::ifstream& ifs) const;

    /// Parse and validate the dimension section
    void parseDimensionSection(std::ifstream& ifs) const;

    /// Parse the nodes section
    void parseNodesSection(std::ifstream& ifs, const Token& token);

    /// Parse a binary nodes data block (single or double precision)
    void parseNodesSectionBinary
    (
        std::ifstream& ifs,
        const Token& token,
        bool doublePrecision
    );

    /// Parse the cells section
    void parseCellsSection(std::ifstream& ifs, const Token& token);

    /// Consume a binary cells data block (cell types are unused)
    void parseCellsSectionBinary(std::ifstream& ifs, const Token& token);

    /// Parse the faces section
    void parseFacesSection(std::ifstream& ifs, const Token& token);

    /// Parse a binary faces data block (int32 connectivity)
    void parseFacesSectionBinary(std::ifstream& ifs, const Token& token);

    /// Parse the boundaries section
    void parseBoundariesSection(std::ifstream& ifs, const Token& token);

    /// Build cell-face connectivity and neighbor relationships
    void buildTopology();

    /// Validate mesh integrity
    void validateMesh() const;
};
