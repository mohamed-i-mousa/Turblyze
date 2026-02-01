/******************************************************************************
 * @file MeshReader.hpp
 * @brief Fluent mesh file reader for ANSYS mesh files
 *
 * This header defines the MeshReader class which reads mesh data from Fluent
 * (.msh) files and converts them into the internal structure
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

#ifndef MESH_READER_HPP
#define MESH_READER_HPP

#include <string>
#include <vector>
#include <fstream>
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
    std::vector<Vector> moveNodes() { return std::move(nodes_); }

    /**
     * @brief Transfer ownership of faces data
     * @return Moved vector of mesh faces
     */
    std::vector<Face> moveFaces() { return std::move(faces_); }

    /**
     * @brief Transfer ownership of cells data
     * @return Moved vector of mesh cells
     */
    std::vector<Cell> moveCells() { return std::move(cells_); }

    /**
     * @brief Transfer ownership of boundary patches data
     * @return Moved vector of boundary patches
     */
    std::vector<BoundaryPatch> moveBoundaryPatches()
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

    /// Comment section "(0"
    static const std::string MSH_COMMENT;
    /// Dimension section "(2"
    static const std::string MSH_DIMENSION;
    /// Nodes section "(10"
    static const std::string MSH_NODES;
    /// Cells section "(12"
    static const std::string MSH_CELLS;
    /// Faces section "(13"
    static const std::string MSH_FACES;
    /// Boundaries section "(45"
    static const std::string MSH_BOUNDARIES;

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
     * 1. Assigns cell indices and clears existing connectivity
     * 2. Maps faces to their owner and neighbor cells
     * 3. Assigns face signs (+1 owner, -1 neighbor)
     * 4. Deduplicates and sorts neighbor cell lists
     */
    void buildTopology();

    /**
     * @brief Validate mesh integrity
     */
    void validateMesh() const;

    /**
     * @brief Print mesh loading summary to stdout
     */
    void printSummary() const;

// Static utility methods

    /**
     * @brief Convert hexadecimal string to size_t
     * @param hexStr Hexadecimal string to convert
     * @return Converted decimal value
     */
    static size_t hexToDec(const std::string& hexStr);

    /**
     * @brief Convert decimal string to size_t
     * @param decStr Decimal string to convert
     * @return Converted decimal value
     */
    static size_t strToDec(const std::string& decStr);

    /**
     * @brief Safely convert 1-based Fluent index to 0-based index
     * @param fluentIdx 1-based index from Fluent file
     * @param context Description of where this index is used (for error messages)
     * @return 0-based index (fluentIdx - 1)
     * @throws std::runtime_error if fluentIdx is 0
     */
    static size_t safeFluentIndexConvert
    (
        size_t fluentIdx,
        const std::string& context
    );
};

#endif // MESH_READER_HPP