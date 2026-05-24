/******************************************************************************
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

#include <string>
#include <string_view>
#include <vector>
#include <iosfwd>

#include "Scalar.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"


class MeshReader
{
public:
    /// Construct MeshReader and parse the given Fluent mesh file
    explicit MeshReader(const std::string& filePath);

// Move accessors

    /// Transfer ownership of nodes data
    [[nodiscard]] std::vector<Vector> moveNodes() noexcept
    {
        return std::move(nodes_);
    }

    /// Transfer ownership of faces data
    [[nodiscard]] std::vector<Face> moveFaces() noexcept
    {
        return std::move(faces_);
    }

    /// Transfer ownership of cells data
    [[nodiscard]] std::vector<Cell> moveCells() noexcept
    {
        return std::move(cells_);
    }

    /// Transfer ownership of boundary patches data
    [[nodiscard]] std::vector<BoundaryPatch> moveBoundaryPatches() noexcept
    {
        return std::move(boundaryPatches_);
    }

private:

// Mesh data

    /// All mesh node coordinates
    std::vector<Vector> nodes_;

    /// All mesh faces
    std::vector<Face> faces_;

    /// All mesh cells
    std::vector<Cell> cells_;

    /// All boundary patches
    std::vector<BoundaryPatch> boundaryPatches_;

// Section identifier constants

    static constexpr std::string_view MSH_COMMENT    = "(0";
    static constexpr std::string_view MSH_DIMENSION  = "(2";
    static constexpr std::string_view MSH_NODES      = "(10";
    static constexpr std::string_view MSH_CELLS      = "(12";
    static constexpr std::string_view MSH_FACES      = "(13";
    static constexpr std::string_view MSH_BOUNDARIES = "(45";

// Private parsing methods

    /// Parse the complete mesh file
    void parseFile(const std::string& filePath);

    /// Parse the comment section and skip its contents
    void parseCommentSection(std::ifstream& ifs) const;

    /// Parse and validate the dimension section
    void parseDimensionSection(std::ifstream& ifs) const;

    /// Parse the nodes section
    void parseNodesSection(std::ifstream& ifs, const std::string& token);

    /// Parse the cells section
    void parseCellsSection(std::ifstream& ifs, const std::string& token);

    /// Parse the faces section
    void parseFacesSection(std::ifstream& ifs, const std::string& token);

    /// Parse the boundaries section
    void parseBoundariesSection(std::ifstream& ifs, const std::string& token);

// Building topology & connectivity methods

    /// Build cell-face connectivity and neighbor relationships
    void buildTopology();

    /// Validate mesh integrity
    void validateMesh() const;

    /// Print mesh loading summary to stdout
    void printSummary() const;

// Static utility methods

    /// Convert hexadecimal string to size_t
    [[nodiscard]] static size_t hexToDec(std::string_view hexStr);

    /// Convert decimal string to size_t
    [[nodiscard]] static size_t strToDec(std::string_view decStr);

    /// Safely convert 1-based Fluent index to 0-based index
    [[nodiscard]] static size_t safeFluentIndexConvert
    (
        size_t fluentIdx,
        std::string_view context
    );

// Fluent BC type mapping

    /// Mapping entry from Fluent type string to enum
    struct BCMapping
    {
        std::string_view fluentType;
        PatchType patchType;
    };

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

    /// Map Fluent boundary type string to enumeration
    [[nodiscard]] static PatchType mapFluentBCToEnum
    (
        std::string_view fluentType
    );
};
