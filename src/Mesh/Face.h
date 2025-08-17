#ifndef FACE_H
#define FACE_H

#include <vector>
#include <string>
#include <optional>

#include "Scalar.h"
#include "Vector.h"    

/**
 * @brief Represents a face in the computational mesh
 * 
 * A face is a boundary between cells or a boundary of the domain.
 * It stores connectivity information, geometric properties, and
 * distance vectors for finite volume calculations.
 */
struct Face 
{
    /// Unique face identifier
    size_t id = 0;
    
    /// Indices of nodes that define this face
    std::vector<size_t> nodeIndices;
    
    /// Index of the owner cell
    size_t ownerCell = 0;
    
    /// Index of neighbor cell (nullopt for boundary faces)
    std::optional<size_t> neighbourCell;

    /// Second moment integrals for centroid calculation
    Scalar x2_integral = 0.0;
    Scalar y2_integral = 0.0;
    Scalar z2_integral = 0.0;

    /// Face geometric center
    Vector centroid;
    
    /// Face normal vector (unit vector)
    Vector normal;
    
    /// Face area
    Scalar area = 0.0;

    /// Distance vector from owner cell center to face center
    Vector d_Pf;
    
    /// Distance vector from neighbor cell center to face center
    std::optional<Vector> d_Nf;
    
    /// Magnitude of d_Pf
    Scalar d_Pf_mag = 0.0;
    
    /// Magnitude of d_Nf
    std::optional<Scalar> d_Nf_mag;
    
    /// Unit vector in d_Pf direction
    Vector e_Pf;
    
    /// Unit vector in d_Nf direction
    std::optional<Vector> e_Nf;

    /// Flag indicating if geometric properties calculated
    bool geometricPropertiesCalculated = false;
    
    /// Flag indicating if distance properties calculated
    bool distancePropertiesCalculated = false;

    /**
     * @brief Default constructor
     */
    Face() = default;

    /**
     * @brief Constructor for internal faces
     * @param faceId Unique face identifier
     * @param nodes Indices of face nodes
     * @param owner Index of owner cell
     * @param neighbour Index of neighbor cell
     */
    Face
    (
        size_t faceId, 
        const std::vector<size_t>& nodes, 
        size_t owner,
        size_t neighbour
    ) : id(faceId), 
        nodeIndices(nodes), 
        ownerCell(owner), 
        neighbourCell(neighbour) {}
    
    /**
     * @brief Constructor for boundary faces
     * @param faceId Unique face identifier
     * @param nodes Indices of face nodes
     * @param owner Index of owner cell
     */
    Face
    (
        size_t faceId, 
        const std::vector<size_t>& nodes, 
        size_t owner
    ) : id(faceId), 
        nodeIndices(nodes), 
        ownerCell(owner),
        neighbourCell(std::nullopt) {}

    /**
     * @brief Calculate geometric properties of the face
     * @param allNodes Vector of all mesh nodes
     * @throws std::out_of_range if node index is invalid
     * @throws std::runtime_error if face is degenerate
     * 
     * Calculates face area, centroid, normal vector, and second moment
     * integrals. For triangles, uses direct cross product. For polygons,
     * decomposes into triangles and uses weighted averaging.
     * Sets geometricPropertiesCalculated flag to true upon success.
     */
    void calculateGeometricProperties(const std::vector<Vector>& allNodes);

    /**
     * @brief Calculate distance properties of the face
     * @tparam CellContainer Container type holding cells
     * @param allCells Container of all mesh cells
     * @throws std::runtime_error if distances are near zero
     * 
     * Calculates distance vectors, magnitudes, and unit vectors
     * from cell centers to face center. For boundary faces,
     * only owner cell distances are calculated.
     */
    template<typename CellContainer>
    void calculateDistanceProperties(const CellContainer& allCells);

    /**
     * @brief Check if this is a boundary face
     * @return True if face is on domain boundary
     */
    bool isBoundary() const
    {
        return !neighbourCell.has_value();
    }
};

/**
 * @brief Stream output operator for Face
 * @param os Output stream
 * @param f Face to output
 * @return Reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Face& f);

#endif