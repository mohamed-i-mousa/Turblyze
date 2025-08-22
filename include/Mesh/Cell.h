#ifndef CELL_H
#define CELL_H

#include <vector>

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"

/**
 * @brief Represents a computational cell in the mesh
 * 
 * A cell is a finite volume element bounded by faces. It stores
 * connectivity information, geometric properties, and provides
 * methods for calculating volume and centroid.
 */
struct Cell
{
    /// Unique cell identifier
    size_t id = 0;
    
    /// Indices of faces that bound this cell
    std::vector<size_t> faceIndices;
    
    /// Indices of neighboring cells
    std::vector<size_t> neighbourCellIndices;
    
    /// Face normal direction signs (+1 outward, -1 inward)
    std::vector<int> faceSigns;
    
    /// Cell geometric center
    Vector centroid;
    
    /// Cell volume
    Scalar volume = 0.0;
    
    /// Flag indicating if geometry has been calculated
    bool geometricPropertiesCalculated = false;

    /**
     * @brief Default constructor
     */
    Cell() = default;

    /**
     * @brief Constructs cell with connectivity data
     * @param cellId Unique cell identifier
     * @param faces Indices of bounding faces
     * @param neighbours Indices of neighboring cells
     * @param signs Face normal direction signs
     */
    Cell
    (
        size_t cellId, 
        const std::vector<size_t>& faces, 
        const std::vector<size_t>& neighbours, 
        const std::vector<int>& signs
    ) : id(cellId),
        faceIndices(faces),
        neighbourCellIndices(neighbours),
        faceSigns(signs) {}

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
};

/**
 * @brief Stream output operator for Cell
 * @param os Output stream
 * @param c Cell to output
 * @return Reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Cell& c);

#endif