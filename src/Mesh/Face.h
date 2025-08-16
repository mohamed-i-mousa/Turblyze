#ifndef FACE_H
#define FACE_H

#include <vector>
#include <string>
#include <optional>

#include "Scalar.h"
#include "Vector.h"    

struct Face {
    // ----- Members ----- //
    size_t id = 0;                                                  
    std::vector<size_t> nodeIndices;                           
    size_t ownerCell = 0;
    std::optional<size_t> neighbourCell;  // No value = boundary face

    Scalar x2_integral = 0.0;
    Scalar y2_integral = 0.0;
    Scalar z2_integral = 0.0;

    Vector centroid;   
    Vector normal;
    Scalar area = 0.0;

    // Pre-computed distance vectors and properties for efficiency
    Vector d_Pf;
    std::optional<Vector> d_Nf;
    Scalar d_Pf_mag = 0.0;
    std::optional<Scalar> d_Nf_mag;
    Vector e_Pf;
    std::optional<Vector> e_Nf;

    bool geometricPropertiesCalculated = false;
    bool distancePropertiesCalculated = false;

    // ----- Constructors ----- //

    Face() = default;

    // Constructor for internal faces
    Face
    (
        size_t faceId, 
        const std::vector<size_t>& nodes, 
        size_t owner,
        size_t neighbour
    ) 
    :   id(faceId), 
        nodeIndices(nodes), 
        ownerCell(owner), 
        neighbourCell(neighbour) 
        {}
    
    // Constructor for boundary faces
    Face
    (
        size_t faceId, 
        const std::vector<size_t>& nodes, 
        size_t owner
    )
    :   id(faceId), 
        nodeIndices(nodes), 
        ownerCell(owner),
        neighbourCell(std::nullopt) 
        {}

    // ----- Member Methods ----- //

    /* Calculate geometric properties of the face
     *
     * Input: allNodes (Vector of all nodes in the mesh and connectivity).
     * Output: necessary geometric properties of the face like area, centroid,
     * normal, etc.
     * 
     * The function calculates the geometric properties of the face based on the
     * number of nodes.
     * 
     * If the face is a triangle, the function calculates the area and normal
     * using the cross product.
     * 
     * If the face is a polygon, the function decomposes the face into triangles
     * and calculates the area and normal using the cross product
     *
     * The function calculates the centroid of the face using the weighted
     * average of the centroids of the triangles
     * 
     * The function calculates the x2_integral, y2_integral, and z2_integral of
     * the face using the weighted average of the x2_integral, y2_integral, and
     * z2_integral of the triangles
     * 
     * The function sets the geometricPropertiesCalculated flag to true and
     * returns it.
     */
    void calculateGeometricProperties(const std::vector<Vector>& allNodes);

    /* Calculate distance properties of the face
     *
     * Input: allCells (Vector of all cells in the mesh).
     * Output: distance vectors, magnitudes, and unit vectors.
     * 
     * The function calculates:
     * - d_Pf: distance vector from owner cell center to face center
     * - d_Nf: distance vector from neighbor cell center to face center 
     * - d_P, d_N: magnitudes of d_Pf and d_Nf
     * - e_Pf, e_Nf: unit vectors of d_Pf and d_Nf
     */
    template<typename CellContainer>
    void calculateDistanceProperties(const CellContainer& allCells);

    bool isBoundary() const
    {
        return !neighbourCell.has_value();
    }
};

// ----- Operator Overloads ----- //

// Forward declaration for operator<<
std::ostream& operator<<(std::ostream& os, const Face& f);

#endif