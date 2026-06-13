/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshReader.h
 * @brief Fluent mesh file reader for ANSYS mesh files
 *
 * @details This header defines the MeshReader class which reads mesh data
 * from Fluent (.msh) files and converts them into the internal structure
 * (Nodes, Faces, Cells). Currently supports 3D unstructured meshes exported
 * from ANSYS Meshing.
 *
 * @class MeshReader
 * - Parsing file sections (Nodes, Cells, Faces, Boundaries)
 * - Constructing the topology and connectivity
 * - Establishing owner-neighbor relationships for all faces
 * - Identifying and grouping boundary patches
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

    /// Mapping entry from Fluent type string to enum
    struct BCMapping
    {
        TokenRef fluentType;
        PatchType patchType;
    };

    /// All mesh node coordinates
    NodeList nodes_;

    /// All mesh faces
    FaceList faces_;

    /// All mesh cells
    CellList cells_;

    /// All boundary patches
    PatchList boundaryPatches_;

    /// Fluent mesh-file section identifier tokens
    static constexpr TokenRef MSH_COMMENT    = "(0";
    static constexpr TokenRef MSH_DIMENSION  = "(2";
    static constexpr TokenRef MSH_NODES      = "(10";
    static constexpr TokenRef MSH_CELLS      = "(12";
    static constexpr TokenRef MSH_FACES      = "(13";
    static constexpr TokenRef MSH_BOUNDARIES = "(45";

    /// Lookup table for Fluent BC type string to enum mapping
    static constexpr BCMapping bcMappings_[] =
    {
        {"velocity-inlet",   PatchType::velocityInlet},
        {"pressure-inlet",   PatchType::pressureInlet},
        {"pressure-outlet",  PatchType::pressureOutlet},
        {"wall",             PatchType::wall},
        {"symmetry",         PatchType::symmetry},
        {"periodic",         PatchType::periodic},
        {"periodic-shadow",  PatchType::periodic},
        {"mass-flow-inlet",  PatchType::massFlowInlet},
        {"outflow",          PatchType::outflow},
        {"interface",        PatchType::interface},
        {"interior",         PatchType::interior},
        {"solid",            PatchType::solid},
        {"fluid",            PatchType::fluid}
    };

// ****************************** Private Methods *****************************

private:

    /// Parse the complete mesh file
    void parseFile(const FilePath& filePath);

    /// Parse the comment section and skip its contents
    void parseCommentSection(std::ifstream& ifs) const;

    /// Parse and validate the dimension section
    void parseDimensionSection(std::ifstream& ifs) const;

    /// Parse the nodes section
    void parseNodesSection(std::ifstream& ifs, TokenRef token);

    /// Parse the cells section
    void parseCellsSection(std::ifstream& ifs, TokenRef token);

    /// Parse the faces section
    void parseFacesSection(std::ifstream& ifs, TokenRef token);

    /// Parse the boundaries section
    void parseBoundariesSection(std::ifstream& ifs, TokenRef token);

    /// Build cell-face connectivity and neighbor relationships
    void buildTopology();

    /// Validate mesh integrity
    void validateMesh() const;

    /// Print mesh loading summary to stdout
    void printSummary() const;

    /// Convert hexadecimal string to Count
    [[nodiscard]] static Count hexToDec(TokenRef hexStr);

    /// Convert decimal string to Count
    [[nodiscard]] static Count strToDec(TokenRef decStr);

    /// Safely convert 1-based Fluent index to 0-based index
    [[nodiscard]] static Index safeFluentIndexConvert
    (
        Count fluentIdx,
        MessageRef context
    );

    /// Map Fluent boundary type string to enumeration
    [[nodiscard]] static PatchType mapFluentBCToEnum
    (
        TokenRef fluentType
    );
};
