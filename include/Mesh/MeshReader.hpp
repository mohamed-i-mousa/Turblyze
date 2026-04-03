/******************************************************************************
 * @file MeshReader.hpp
 * @brief Fluent mesh file reader for ANSYS mesh files
 *
 * @details This header defines the MeshReader class which reads mesh data
 * from Fluent (.msh) files and converts them into the internal structure
 * (Nodes, Faces, Cells). Currently supports 3D unstructured meshes exported
 * from ANSYS Meshing.
 *
 * @class MeshReader
 *
 * The mesh reading process involves:
 * - Parsing file sections (Nodes, Cells, Faces, Boundaries)
 * - Constructing the topology and connectivity
 * - Establishing owner-neighbor relationships for all faces
 * - Identifying and grouping boundary patches
 *
 * Supported Fluent face types (hexadecimal):
 * - "2" = internal, "3" = wall, "4" = pressure-inlet
 * - "5" = pressure-outlet, "7" = symmetry, "8" = periodic-shadow
 * - "9" = pressure-far-field, "a" = velocity-inlet, "c" = periodic
 * - "e" = fan/porous-jump, "14" = mass-flow-inlet, "18" = interface
 * - "1F" = parent, "24" = outflow, "25" = axis
 *
 * Supported Fluent element types (hexadecimal):
 * - "0" = mixed, "2" = line/edge, "3" = triangular
 * - "4" = quadrilateral, "5" = polygonal
 *****************************************************************************/

#pragma once

#include <string>
#include <string_view>
#include <vector>
#include <iosfwd>

#include "Scalar.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "BoundaryPatch.hpp"


class MeshReader
{
public:
    /**
     * @brief Construct MeshReader and parse the given Fluent mesh file
     * @param filePath Path to the Fluent mesh file (.msh format)
     */
    explicit MeshReader(const std::string& filePath);

// Move accessors

    /**
     * @brief Transfer ownership of nodes data
     * @return Moved vector of node coordinates
     */
    [[nodiscard("Moved mesh data must be captured")]]
    std::vector<Vector> moveNodes() noexcept
    {
        return std::move(nodes_);
    }

    /**
     * @brief Transfer ownership of faces data
     * @return Moved vector of mesh faces
     */
    [[nodiscard("Moved mesh data must be captured")]]
    std::vector<Face> moveFaces() noexcept
    {
        return std::move(faces_);
    }

    /**
     * @brief Transfer ownership of cells data
     * @return Moved vector of mesh cells
     */
    [[nodiscard("Moved mesh data must be captured")]]
    std::vector<Cell> moveCells() noexcept
    {
        return std::move(cells_);
    }

    /**
     * @brief Transfer ownership of boundary patches data
     * @return Moved vector of boundary patches
     */
    [[nodiscard("Moved mesh data must be captured")]]
    std::vector<BoundaryPatch> moveBoundaryPatches() noexcept
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

    /**
     * @brief Parse the complete mesh file
     * @param filePath Path to the Fluent mesh file
     */
    void parseFile(const std::string& filePath);

    /**
     * @brief Parse the comment section and skip its contents
     * @param ifs Input file stream positioned at the comment section
     */
    void parseCommentSection(std::ifstream& ifs) const;

    /**
     * @brief Parse and validate the dimension section
     * @param ifs Input file stream positioned at the dimension section
     */
    void parseDimensionSection(std::ifstream& ifs) const;

    /**
     * @brief Parse the nodes section (header or data block)
     * @param ifs Input file stream positioned after the section identifier
     * @param token The current token being parsed
     */
    void parseNodesSection(std::ifstream& ifs, const std::string& token);

    /**
     * @brief Parse the cells section (header only, allocates storage)
     * @param ifs Input file stream positioned after the section identifier
     * @param token The current token being parsed
     */
    void parseCellsSection(std::ifstream& ifs, const std::string& token);

    /**
     * @brief Parse the faces section (header or connectivity data)
     * @param ifs Input file stream positioned after the section identifier
     * @param token The current token being parsed
     */
    void parseFacesSection(std::ifstream& ifs, const std::string& token);

    /**
     * @brief Parse the boundaries section (zone names and types)
     * @param ifs Input file stream positioned after the section identifier
     * @param token The current token being parsed
     */
    void parseBoundariesSection(std::ifstream& ifs, const std::string& token);

// Building topology & connectivity methods

    /**
     * @brief Build cell-face connectivity and neighbor relationships
     *
     * @details
     * 1. Assigns cell indices and clears existing connectivity
     * 2. Maps faces to their owner and neighbor cells
     * 3. Assigns face signs (+1 owner, -1 neighbor)
     * 4. Deduplicates and sorts neighbor cell lists
     */
    void buildTopology();

    /// Validate mesh integrity
    void validateMesh() const;

    /// Print mesh loading summary to stdout
    void printSummary() const;

// Static utility methods

    /**
     * @brief Convert hexadecimal string to size_t
     * @param hexStr Hexadecimal string to convert
     * @return Converted decimal value
     */
    static size_t hexToDec(std::string_view hexStr);

    /**
     * @brief Convert decimal string to size_t
     * @param decStr Decimal string to convert
     * @return Converted decimal value
     */
    static size_t strToDec(std::string_view decStr);

    /**
     * @brief Safely convert 1-based Fluent index to 0-based index
     * @param fluentIdx 1-based index from Fluent file
     * @param context Where this index is used (in the Error message)
     * @return 0-based index (fluentIdx - 1)
     * @throws std::runtime_error if fluentIdx is 0
     */
    static size_t safeFluentIndexConvert
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
        {"velocity-inlet",   PatchType::VELOCITY_INLET},
        {"pressure-inlet",   PatchType::PRESSURE_INLET},
        {"pressure-outlet",  PatchType::PRESSURE_OUTLET},
        {"wall",             PatchType::WALL},
        {"symmetry",         PatchType::SYMMETRY},
        {"periodic",         PatchType::PERIODIC},
        {"periodic-shadow",  PatchType::PERIODIC},
        {"mass-flow-inlet",  PatchType::MASS_FLOW_INLET},
        {"outflow",          PatchType::OUTFLOW},
        {"interface",        PatchType::INTERFACE},
        {"interior",         PatchType::INTERIOR},
        {"solid",            PatchType::SOLID},
        {"fluid",            PatchType::FLUID}
    };

    /**
     * @brief Maps Fluent boundary type string to enumeration
     * @param fluentType String representation from Fluent mesh file
     * @return Corresponding PatchType enumeration
     */
    [[nodiscard("Computed enum mapping is required")]]
    static PatchType mapFluentBCToEnum
    (
        std::string_view fluentType
    );
};
