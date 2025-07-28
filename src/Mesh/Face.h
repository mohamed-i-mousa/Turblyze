#ifndef FACE_H
#define FACE_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>      // for std::setprecision, std::fixed in operator<<
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

    bool geometricPropertiesCalculated = false;

    // Default constructor >>> Constructor for internal faces >>> Constructor for boundary faces
    Face() = default;

    Face(size_t faceId, const std::vector<size_t>& nodes, size_t owner, size_t neighbour)
        : id(faceId), nodeIndices(nodes), ownerCell(owner), neighbourCell(neighbour) {}
    
    Face(size_t faceId, const std::vector<size_t>& nodes, size_t owner)
        : id(faceId), nodeIndices(nodes), ownerCell(owner), neighbourCell(std::nullopt) {}

    // ----- Member Methods ----- //

    /* Calculate geometric properties of the face
     * Input: allNodes >>> Vector of all nodes in the mesh
     * Output: necessary geometric properties of the face like area, centroid, normal, etc.
     * The function starts with basic validation of the face's nodes count and node indices
     * The function then calculates the geometric properties of the face based on the number of nodes
     * If the face is a triangle, the function calculates the area and normal using the cross product
     * If the face is a polygon, the function decomposes the face into triangles and calculates the area and normal using the cross product
     * The function then calculates the centroid of the face using the weighted average of the centroids of the triangles
     * The function then calculates the x2_integral, y2_integral, and z2_integral of the face using the weighted average of the x2_integral, y2_integral, and z2_integral of the triangles
     * The function then sets the geometricPropertiesCalculated flag to true
     * The function then returns the geometricPropertiesCalculated flag
     */
    void calculateGeometricProperties(const std::vector<Vector>& allNodes) {
        geometricPropertiesCalculated = false;         
        const size_t nNodes = nodeIndices.size();

        if (nNodes < 3) {
            throw std::runtime_error("Warning: Face " + std::to_string(id) + " has fewer than 3 nodes (" + std::to_string(nNodes)
                      + "). Cannot calculate geometric properties. Skipping.");
        }

        for (size_t i = 0; i < nNodes; ++i) {
            if (nodeIndices[i] >= allNodes.size()) {
                throw std::out_of_range("Error calculating properties for Face " + std::to_string(id) +
                                        ": Node index " + std::to_string(nodeIndices[i]) + " out of range (" +
                                        "Node list size: " + std::to_string(allNodes.size()) + ").");
            }
        }

        // CASE 1: Face is "Triangle" (nNodes == 3)
        if (nNodes == 3) {
            const Vector& p1 = allNodes[nodeIndices[0]];
            const Vector& p2 = allNodes[nodeIndices[1]];
            const Vector& p3 = allNodes[nodeIndices[2]];

            centroid = (p1 + p2 + p3) / S(3.0);

            Vector vecA = p2 - p1;
            Vector vecB = p3 - p1;

            Vector crossProd = cross(vecB, vecA);
            Scalar crossProdMag = crossProd.magnitude();

            if (std::abs(crossProdMag) < AREA_TOLERANCE) {
                throw std::runtime_error("Face " + std::to_string(id) + " is degenerate.");
            }

            else {
                area = S(0.5) * crossProdMag;
                if (std::abs(crossProdMag) > DIVISION_TOLERANCE) {
                    normal = crossProd / crossProdMag;
                } else {
                    normal = Vector(S(0.0), S(0.0), S(0.0));
                }

                x2_integral = (p1.x*p1.x + p2.x*p2.x + p3.x*p3.x + p1.x*p2.x + p1.x*p3.x + p2.x*p3.x) * area / S(6.0);
                y2_integral = (p1.y*p1.y + p2.y*p2.y + p3.y*p3.y + p1.y*p2.y + p1.y*p3.y + p2.y*p3.y) * area / S(6.0);
                z2_integral = (p1.z*p1.z + p2.z*p2.z + p3.z*p3.z + p1.z*p2.z + p1.z*p3.z + p2.z*p3.z) * area / S(6.0);

                geometricPropertiesCalculated = true;
            }
        }
        // CASE 2: Face is "Polygon" (nNodes > 3)
        else {
            x2_integral = 0.0;
            y2_integral = 0.0;
            z2_integral = 0.0;

            Vector faceCenter(0.0, 0.0, 0.0);

            for (size_t i = 0; i < nNodes; ++i) {
                faceCenter += allNodes[nodeIndices[i]];
            }
            if (nNodes > 0) {
                faceCenter /= S(nNodes);
            }

            Scalar totalArea = 0.0;
            Vector weightedCentroidSum(0.0, 0.0, 0.0);
            Vector normalSum(0.0, 0.0, 0.0);

            for (size_t i = 0; i < nNodes; ++i) {
                const Vector& p_i = allNodes[nodeIndices[i]];
                const Vector& p_next = allNodes[nodeIndices[(i + 1) % nNodes]]; // Handles the end point

                const Vector& p1_tri = faceCenter;
                const Vector& p2_tri = p_i;
                const Vector& p3_tri = p_next;

                Vector vecA_tri = p2_tri - p1_tri;
                Vector vecB_tri = p3_tri - p1_tri;

                Vector crossProd_tri = cross(vecB_tri, vecA_tri);
                Scalar triangleArea = S(0.5) * crossProd_tri.magnitude();
                normalSum += crossProd_tri;

                Scalar x2_part = (p1_tri.x * p1_tri.x + p2_tri.x * p2_tri.x + p3_tri.x * p3_tri.x +
                                p1_tri.x * p2_tri.x + p1_tri.x * p3_tri.x + p2_tri.x * p3_tri.x) *
                               triangleArea / S(6.0);
                Scalar y2_part = (p1_tri.y * p1_tri.y + p2_tri.y * p2_tri.y + p3_tri.y * p3_tri.y +
                                p1_tri.y * p2_tri.y + p1_tri.y * p3_tri.y + p2_tri.y * p3_tri.y) *
                               triangleArea / S(6.0);
                Scalar z2_part = (p1_tri.z * p1_tri.z + p2_tri.z * p2_tri.z + p3_tri.z * p3_tri.z +
                                p1_tri.z * p2_tri.z + p1_tri.z * p3_tri.z + p2_tri.z * p3_tri.z) *
                               triangleArea / S(6.0);

                x2_integral += x2_part;
                y2_integral += y2_part;
                z2_integral += z2_part;

                if (triangleArea > AREA_TOLERANCE) {
                    Vector triangleCentroid = (p1_tri + p2_tri + p3_tri) / S(3.0);
                    totalArea += triangleArea;
                    weightedCentroidSum += triangleCentroid * triangleArea;
                }
            }

            area = totalArea;

            if (area < AREA_TOLERANCE) {
                 throw std::runtime_error("Warning: Polygonal Face " + std::to_string(id) + " has near-zero total area. Setting area=0, centroid/normal=(0,0,0).");
            } else {
                 if (std::abs(area) > DIVISION_TOLERANCE) {
                     centroid = weightedCentroidSum / area;
                 } else {
                     throw std::runtime_error("Warning: Polygonal Face " + std::to_string(id) + " has near-zero total area.");
                 }
                 normal = normalSum.normalized();
                 geometricPropertiesCalculated = true; 
            }
        }
    }


    bool isBoundary() const {
        return !neighbourCell.has_value();
    }
};

// ----- Operator Overloads (Non-Member Methods) ----- //

/* This function is used to print the face to the console
 * It prints the face's id, nodes, owner, neighbour, and calculated properties if available
 * The function sets the precision for floating Vector output within this scope
 * The function then prints the face's id, nodes, owner, neighbour, and calculated properties if available
 * The function then restores the precision for floating Vector output within this scope
 * The function then returns the ostream object
 */
 inline std::ostream& operator<<(std::ostream& os, const Face& f) {
    os << "Face(ID: " << f.id
    << ", Nodes: [";
    for (size_t i = 0; i < f.nodeIndices.size(); ++i) {
        os << f.nodeIndices[i] << (i == f.nodeIndices.size() - 1 ? "" : ", ");
    }
    os << "], Owner: " << f.ownerCell
       << ", Neighbour: " << (f.isBoundary() ? "Boundary" : std::to_string(f.neighbourCell.value_or(0))); // Display "boundary" for boundary

   // Print calculated properties if available
    if (f.geometricPropertiesCalculated) {
        // Set precision for floating Vector output within this scope
        std::ios_base::fmtflags flags = os.flags();
        int                      prec = os.precision();
        os << std::fixed << std::setprecision(6);
        os << ", Centroid: " << f.centroid
        << ", Area: " << f.area
        << ", Normal: " << f.normal;
        os.flags(flags);
        os.precision(prec);
    } else {
        os << ", Geometry: N/A";
    }
    os << ")";
    return os;
}

#endif
