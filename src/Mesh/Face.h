#ifndef FACE_H
#define FACE_H

#include <vector>     // for std::vector
#include <stdexcept>  // for std::out_of_range and std::runtime_error
#include <cmath>      // for std::abs
#include <iostream>   // for std::ostream, std::cerr
#include <iomanip>    // for std::setprecision, std::std::fixed in operator<<
#include <string>     // for std::to_string in operator<<
#include <optional>   // for std::optional

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
    Vector normal;              // Normal vector scaled by face area (Vectors owner->neighbor)
    Scalar area = 0.0;          // Magnitude of the area vector

    bool geometricPropertiesCalculated = false;     // Flag to check if calculations were successful

    // ----- Constructors ----- //

    // Default constructor (useful for creating vectors of Faces before knowing details)
    Face() = default;

    // Constructor for internal faces
    Face(size_t faceId, const std::vector<size_t>& nodes, size_t owner, size_t neighbour)
        : id(faceId), nodeIndices(nodes), ownerCell(owner), neighbourCell(neighbour) {}
    
    // Constructor for boundary faces
    Face(size_t faceId, const std::vector<size_t>& nodes, size_t owner)
        : id(faceId), nodeIndices(nodes), ownerCell(owner), neighbourCell(std::nullopt) {}

    // ----- Member Methods ----- //

    void calculateGeometricProperties(const std::vector<Vector>& allNodes) {
        geometricPropertiesCalculated = false;         
        const size_t nNodes = nodeIndices.size();

        // Basic validation
        if (nNodes < 3) {
            throw std::runtime_error("Warning: Face " + std::to_string(id) + " has fewer than 3 nodes (" + std::to_string(nNodes)
                      + "). Cannot calculate geometric properties. Skipping.");
        }

        // Node index validation
        for (size_t i = 0; i < nNodes; ++i) {
            if (nodeIndices[i] >= allNodes.size()) {
                throw std::out_of_range("Error calculating properties for Face " + std::to_string(id) +
                                        ": Node index " + std::to_string(nodeIndices[i]) + " out of range (" +
                                        "Node list size: " + std::to_string(allNodes.size()) + ").");
            }
        }

        // Create a local vector of points for easier access
        std::vector<Vector> points(nNodes);
        for (size_t i = 0; i < nNodes; ++i) {
            points[i] = allNodes[nodeIndices[i]];
        }

        // CASE 1: Face is "Triangle" (nNodes == 3)
        if (nNodes == 3) {
            const Vector& p1 = points[0];
            const Vector& p2 = points[1];
            const Vector& p3 = points[2];

            // Calculate the Centroid
            centroid = (p1 + p2 + p3) / S(3.0);

            // Calculate area & normal (using cross product)
            // Use consistent ordering to ensure normal direction
            Vector vecA = p2 - p1;
            Vector vecB = p3 - p1;
            Vector crossProd = cross(vecB, vecA);
            Scalar crossProdMag = crossProd.magnitude();

            // Check the minimum area limit 
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

            // Calculate face center (vertices average)
            Vector faceCenter(0.0, 0.0, 0.0);
            for (const auto& p : points) {
                faceCenter += p;
            }
            if (nNodes > 0) {
                faceCenter /= S(nNodes);
            } 
            // Decompose into triangles (Center, P_i, P_{i+1}) and sum properties
            Scalar totalArea = 0.0;
            Vector weightedCentroidSum(0.0, 0.0, 0.0);
            Vector normalSum(0.0, 0.0, 0.0);

            for (size_t i = 0; i < nNodes; ++i) {
                const Vector& p_i = points[i];
                const Vector& p_next = points[(i + 1) % nNodes]; // Handles the end point

                // Calculate area and centroid for triangles 
                const Vector& p1_tri = faceCenter;
                const Vector& p2_tri = p_i;
                const Vector& p3_tri = p_next;

                Vector vecA_tri = p2_tri - p1_tri;
                Vector vecB_tri = p3_tri - p1_tri;

                Vector crossProd_tri = cross(vecB_tri, vecA_tri);
                Scalar triangleArea = S(0.5) * crossProd_tri.magnitude();
                normalSum += crossProd_tri;          // inside the loop

                Scalar x2_part = (p1_tri.x * p1_tri.x + p2_tri.x * p2_tri.x + p3_tri.x * p3_tri.x +
                                p1_tri.x * p2_tri.x + p1_tri.x * p3_tri.x + p2_tri.x * p3_tri.x) *
                               triangleArea / S(6.0);
                Scalar y2_part = (p1_tri.y * p1_tri.y + p2_tri.y * p2_tri.y + p3_tri.y * p3_tri.y +
                                p1_tri.y * p2_tri.y + p1_tri.y * p3_tri.y + p2_tri.y * p3_tri.y) *
                               triangleArea / S(6.0);
                Scalar z2_part = (p1_tri.z * p1_tri.z + p2_tri.z * p2_tri.z + p3_tri.z * p3_tri.z +
                                p1_tri.z * p2_tri.z + p1_tri.z * p3_tri.z + p2_tri.z * p3_tri.z) *
                               triangleArea / S(6.0);

                // Accumulate into face integrals
                x2_integral += x2_part;
                y2_integral += y2_part;
                z2_integral += z2_part;

                // Accumulate area and weighted centroid
                if (triangleArea > AREA_TOLERANCE) {
                    Vector triangleCentroid = (p1_tri + p2_tri + p3_tri) * S(1.0 / 3.0);
                    totalArea += triangleArea;
                    weightedCentroidSum += triangleCentroid * triangleArea;
                }
            }

            // Final polygon area
            area = totalArea;

            // Final polygon centroid
            if (area < AREA_TOLERANCE) {
                 throw std::runtime_error("Warning: Polygonal Face " + std::to_string(id) + " has near-zero total area. Setting area=0, centroid/normal=(0,0,0).");

            } else {
                 if (std::abs(area) > DIVISION_TOLERANCE) {
                     centroid = weightedCentroidSum / area;
                 } else {
                     centroid = Vector(S(0.0), S(0.0), S(0.0));
                 }
                 normal = normalSum.normalized();     // after the loop
                 geometricPropertiesCalculated = true; 
            }
        }
    }

    bool isBoundary() const {
        return !neighbourCell.has_value();
    }

    size_t getNeighbour() const {
        if (!neighbourCell.has_value()) {
            throw std::runtime_error("Face " + std::to_string(id) + " is a boundary face");
        }
        return neighbourCell.value();
    }
};

// ----- Operator Overloads (Non-Member Methods) ----- //

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
        std::ios_base::fmtflags flags = os.flags();         // Save current flags
        int                      prec = os.precision();     // Save current precision
        os << std::fixed << std::setprecision(6);           // Set precision 
        os << ", Centroid: " << f.centroid
        << ", Area: " << f.area
        << ", Normal: " << f.normal;
        os.flags(flags);                                    // Restore flags
        os.precision(prec);                                 // Restore precision
    } else {
        os << ", Geometry: N/A";
    }
    os << ")";
    return os;
}

#endif
