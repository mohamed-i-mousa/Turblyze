#ifndef CELL_H
#define CELL_H

#include <vector>

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"


struct Cell
{
    // ----- Members ----- //
    size_t id = 0;
    std::vector<size_t> faceIndices;
    std::vector<size_t> neighbourCellIndices;
    std::vector<int> faceSigns;                   // to adjust normal vector direction to be always pointing outward 
    Vector centroid;
    Scalar volume = 0.0;
    bool geometricPropertiesCalculated = false;

    // ----- Constructors ----- //

    Cell() = default;

    Cell
    (
        size_t cellId, 
        const std::vector<size_t>& faces, 
        const std::vector<size_t>& neighbours, 
        const std::vector<int>& signs
    )
    : id(cellId),
      faceIndices(faces),
      neighbourCellIndices(neighbours),
      faceSigns(signs)
      {}

    /* Calculate geometric properties of the cell
     *
     * Input: allFaces.
     * Output: Cell volume and centroid.
     * 
     * The function calculates the geometric properties of the cell based on
     * the number of faces.
     * The cell volume is calculated using the divergence theorem: 
     *      V = (1/3) * sum(face_centroid . face_area_vector)
     * 
     * The cell centroid is calculated using the second moments of the faces.
     *
     * The function sets the geometricPropertiesCalculated flag to true and
     * returns it.
     */
    void calculateGeometricProperties(const std::vector<Face>& allFaces);
};

// Forward declaration for operator<<
std::ostream& operator<<(std::ostream& os, const Cell& c);

#endif