/******************************************************************************
 * @file Cell.h
 * @brief Represents a computational cell in the mesh
 * 
 * @class Cell 
 * 
 * A cell is a finite volume element bounded by faces. It stores
 * connectivity information, geometric properties, and provides
 * methods for calculating volume and centroid.
 *****************************************************************************/

#ifndef CELL_H
#define CELL_H

#include <vector>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"


class Cell
{
public:

    /// Default constructor 
    Cell() = default;

    /**
     * @brief Constructs cell with connectivity data
     * @param cellId Unique cell identifier
     * @param faces Indices of bounding faces
     * @param neighbors Indices of neighboring cells
     * @param signs Face normal direction signs
     */
    Cell
    (
        size_t cellId, 
        const std::vector<size_t>& faces, 
        const std::vector<size_t>& neighbors, 
        const std::vector<int>& signs
    ) : id_(cellId),
        faceIndices_(faces),
        neighborCellIndices_(neighbors),
        faceSigns_(signs) {}

// Setter methods
    
    /** 
     * @brief Set cell identifier
     * @param cellId Unique cell ID
     */
    void setId(size_t cellId) { id_ = cellId; }
    
    /** 
     * @brief Add face index to cell connectivity
     * @param faceIdx Index of bounding face
     */
    void addFaceIndex(size_t faceIdx) { faceIndices_.push_back(faceIdx); }
    
    /** 
     * @brief Add face normal direction sign
     * @param sign Direction sign (+1 outward, -1 inward)
     */
    void addFaceSign(int sign) { faceSigns_.push_back(sign); }
    
    /** 
     * @brief Clear all face indices
     */
    void clearFaceIndices() { faceIndices_.clear(); }
    
    /** 
     * @brief Add neighboring cell index
     * @param neighborIdx Index of neighboring cell
     */
    void addNeighborCellIndex(size_t neighborIdx)
    {
        neighborCellIndices_.push_back(neighborIdx);
    }
    
    /** 
     * @brief Set all neighbor cell indices
     * @param neighbors Vector of neighboring cell indices
     */
    void setNeighborCellIndices(const std::vector<size_t>& neighbors)
    {
        neighborCellIndices_ = neighbors;
    }
    
    /** 
     * @brief Clear all neighbor cell indices
     */
    void clearNeighborCellIndices() { neighborCellIndices_.clear(); }

// Accessor methods

    /** 
     * @brief Get cell identifier 
     * @return Unique cell ID 
     */
    size_t id() const { return id_; }

    /** 
     * @brief Get bounding face indices 
     * @return Vector of face indices 
     */
    const std::vector<size_t>& faceIndices() const { return faceIndices_; }
    
    /** 
     * @brief Get neighboring cell indices 
     * @return Vector of neighbor cell indices 
     */
    const std::vector<size_t>& neighborCellIndices() const
    {
        return neighborCellIndices_;
    }
    
    /** 
     * @brief Get face normal direction signs 
     * @return Vector of signs (+1/-1) 
     */
    const std::vector<int>& faceSigns() const { return faceSigns_; }

    /** 
     * @brief Get cell centroid 
     * @return Cell center coordinates 
     */
    const Vector& centroid() const { return centroid_; }
    
    /** 
     * @brief Get cell volume 
     * @return Cell volume magnitude 
     */
    Scalar volume() const { return volume_; }

    /** 
     * @brief Check if geometric properties calculated 
     * @return True if geometry computed 
     */
    bool geometricPropertiesCalculated() const 
    {
        return geometricPropertiesCalculated_;
    }

    /**
     * @brief Calculate geometric properties of the cell
     * @param allFaces Vector containing all mesh faces
     * 
     * Calculates cell volume using the divergence theorem:
     * V = (1/3) * Σ(face_centroid · face_area_vector)
     * 
     * Calculates cell centroid using second moments of the faces.
     * Sets geometricPropertiesCalculated flag to true upon success.
     */
    void calculateGeometricProperties(const std::vector<Face>& allFaces);

private:
    /// Unique cell identifier
    size_t id_ = 0;
    
    /// Indices of faces that bound this cell
    std::vector<size_t> faceIndices_;
    
    /// Indices of neighboring cells
    std::vector<size_t> neighborCellIndices_;
    
    /// Face normal direction signs (+1 outward, -1 inward)
    std::vector<int> faceSigns_;
    
    /// Cell geometric center
    Vector centroid_;
    
    /// Cell volume
    Scalar volume_ = 0.0;
    
    /// Flag indicating if geometry has been calculated
    bool geometricPropertiesCalculated_ = false;

    friend std::ostream& operator<<(std::ostream& os, const Cell& c);
};

/**
 * @brief Stream output operator for Cell
 * @param os Output stream
 * @param c Cell to output
 * @return Reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Cell& c);

#endif