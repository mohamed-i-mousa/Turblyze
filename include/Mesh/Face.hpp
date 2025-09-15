/******************************************************************************
 * @file Face.h 
 * @brief Represents a face in the computational mesh
 * 
 * @class Face
 * 
 * A face is a boundary between cells or a boundary of the domain.
 * It stores connectivity information, geometric properties, and
 * distance vectors for finite volume calculations.
 *****************************************************************************/

#ifndef FACE_H
#define FACE_H

#include <vector>
#include <string>
#include <optional>

#include "Scalar.hpp"
#include "Vector.hpp"    

class Face 
{
public:

    /// Default constructor 
    Face() = default;

    /**
     * @brief Constructor for internal faces
     * @param faceId Unique face identifier
     * @param nodes Indices of face nodes
     * @param owner Index of owner cell
     * @param neighbor Index of neighbor cell
     */
    Face
    (
        size_t faceId, 
        const std::vector<size_t>& nodes, 
        size_t owner,
        size_t neighbor
    ) : id_(faceId), 
        nodeIndices_(nodes), 
        ownerCell_(owner), 
        neighborCell_(neighbor) {}
    
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
    ) : id_(faceId), 
        nodeIndices_(nodes), 
        ownerCell_(owner),
        neighborCell_(std::nullopt) {}

// Setter methods
    
    /** 
     * @brief Set face identifier
     * @param faceId Unique face ID
     */
    void setId(size_t faceId) { id_ = faceId; }
    
    /** 
     * @brief Set owner cell index
     * @param owner Index of owner cell
     */
    void setOwnerCell(size_t owner) { ownerCell_ = owner; }
    
    /** 
     * @brief Set neighbor cell index
     * @param neighbor Index of neighbor cell
     */
    void setNeighborCell(size_t neighbor) { neighborCell_ = neighbor; }
    
    /** 
     * @brief Set neighbor cell to null (boundary face)
     */
    void setNeighborCell(std::nullopt_t) { neighborCell_ = std::nullopt; }
    
    /** 
     * @brief Add node index to face connectivity
     * @param nodeIdx Index of node to add
     */
    void addNodeIndex(size_t nodeIdx) { nodeIndices_.push_back(nodeIdx); }
    
    /** 
     * @brief Clear all node indices
     */
    void clearNodeIndices() { nodeIndices_.clear(); }

// Accessor methods

    /** 
     * @brief Get face identifier 
     * @return Unique face ID 
     */
    size_t id() const { return id_; }
    
    /** 
     * @brief Get node connectivity 
     * @return Vector of node indices 
     */
    const std::vector<size_t>& nodeIndices() const { return nodeIndices_; }
    
    /** 
     * @brief Get owner cell index 
     * @return Index of owner cell 
     */
    size_t ownerCell() const { return ownerCell_; }
    
    /** 
     * @brief Get neighbor cell index 
     * @return Optional neighbor cell index 
     */
    const std::optional<size_t>& neighborCell() const { return neighborCell_; }

    /** 
     * @brief Get face centroid 
     * @return Face center coordinates 
     */
    const Vector& centroid() const { return centroid_; }
    
    /** 
     * @brief Get face normal vector 
     * @return Unit normal vector 
     */
    const Vector& normal() const { return normal_; }
    
    /** 
     * @brief Get face area 
     * @return Face area magnitude 
     */
    Scalar area() const { return area_; }
    
    /** 
     * @brief Get second moment integrals for centroid calculation 
     */
    Scalar x2_integral() const { return x2_integral_; }
    Scalar y2_integral() const { return y2_integral_; }
    Scalar z2_integral() const { return z2_integral_; }

    /** 
     * @brief Get owner cell distance vector 
     * @return Vector from owner to face 
     */
    const Vector& d_Pf() const { return d_Pf_; }
    
    /** 
     * @brief Get neighbor cell distance vector 
     * @return Optional vector from neighbor to face
     */
    const std::optional<Vector>& d_Nf() const { return d_Nf_; }
    
    /** 
     * @brief Get owner cell distance magnitude 
     * @return Distance from owner to face 
     */
    Scalar d_Pf_mag() const { return d_Pf_mag_; }
    
    /** 
     * @brief Get neighbor cell distance magnitude 
     * @return Optional distance from neighbor to face
     */
    const std::optional<Scalar>& d_Nf_mag() const { return d_Nf_mag_; }
    
    /** 
     * @brief Get owner cell unit vector 
     * @return Unit vector from owner to face
     */
    const Vector& e_Pf() const { return e_Pf_; }
    
    /** 
     * @brief Get neighbor cell unit vector 
     * @return Optional unit vector from neighbor to face 
     */
    const std::optional<Vector>& e_Nf() const { return e_Nf_; }

    /** 
     * @brief Check if geometric properties calculated 
     * @return True if geometry computed 
     */
    bool geometricPropertiesCalculated() const 
    {
        return geometricPropertiesCalculated_; 
    }
    
    /** 
     * @brief Check if distance properties calculated 
     * @return True if distances computed 
     */
    bool distancePropertiesCalculated() const
    {
        return distancePropertiesCalculated_;
    }

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
        return !neighborCell_.has_value();
    }

private:
    /// Unique face identifier
    size_t id_ = 0;
    
    /// Indices of nodes that define this face
    std::vector<size_t> nodeIndices_;
    
    /// Index of the owner cell
    size_t ownerCell_ = 0;
    
    /// Index of neighbor cell (nullopt for boundary faces)
    std::optional<size_t> neighborCell_;

    /// Second moment integrals for centroid calculation
    Scalar x2_integral_ = 0.0;
    Scalar y2_integral_ = 0.0;
    Scalar z2_integral_ = 0.0;

    /// Face geometric center
    Vector centroid_;
    
    /// Face normal vector (unit vector)
    Vector normal_;
    
    /// Face area
    Scalar area_ = 0.0;

    /// Distance vector from owner cell center to face center
    Vector d_Pf_;
    
    /// Distance vector from neighbor cell center to face center
    std::optional<Vector> d_Nf_;
    
    /// Magnitude of d_Pf
    Scalar d_Pf_mag_ = 0.0;
    
    /// Magnitude of d_Nf
    std::optional<Scalar> d_Nf_mag_;
    
    /// Unit vector in d_Pf direction
    Vector e_Pf_;
    
    /// Unit vector in d_Nf direction
    std::optional<Vector> e_Nf_;

    /// Flag indicating if geometric properties calculated
    bool geometricPropertiesCalculated_ = false;
    
    /// Flag indicating if distance properties calculated
    bool distancePropertiesCalculated_ = false;

    /// friend function for operator << 
    friend std::ostream& operator<<(std::ostream& os, const Face& f);
};

std::ostream& operator<<(std::ostream& os, const Face& f);

#endif