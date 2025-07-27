#ifndef CELL_H
#define CELL_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"


struct Cell {
    // ----- Members ----- //
    size_t id = 0;
    std::vector<size_t> faceIndices;
    std::vector<size_t> neighbourCellIndices;
    std::vector<int> faceSigns;                   // This to adjust normal vector direction to be always pointing outward 
    Vector centroid;
    Scalar volume = 0.0;
    bool geometricPropertiesCalculated = false;

    // Default constructor >>> Constructor for internal faces >>> Constructor for boundary faces

    Cell() = default;

    Cell(size_t cellId, const std::vector<size_t>& faces, const std::vector<size_t>& neighbours, const std::vector<int>& signs)
        : id(cellId),
          faceIndices(faces),
          faceSigns(signs),
          neighbourCellIndices(neighbours)
          {}

    /* Calculate geometric properties of the cell
     * Input: allFaces >>> Vector of all faces in the mesh
     * Output: necessary geometric properties of the cell like volume, centroid, etc.
     * The function starts with basic validation of the cell's faces count and face indices
     * The function then calculates the geometric properties of the cell based on the number of faces
     * The cell volume is calculated using the divergence theorem: V = (1/3) * sum(face_centroid . face_area_vector)
     * The cell centroid is calculated using the second moments of the faces
     */
    void calculateGeometricProperties(const std::vector<Face>& allFaces) {
        geometricPropertiesCalculated = false;
        volume = 0.0;
        centroid = Vector(0.0, 0.0, 0.0);
        Vector centroidSum(0.0, 0.0, 0.0);

        if (faceIndices.empty()) {
            throw std::runtime_error("Warning: Cell " + std::to_string(id) + " has no assigned faces. Calculations cannot proceed!");
        }

        for (size_t i = 0; i < faceIndices.size(); ++i) {
            size_t faceIndex = faceIndices[i];

            // Index validation
            if (faceIndex >= allFaces.size()) {
                throw std::out_of_range("Error in Cell" + std::to_string(id) + "calculation: " +
                "Face index " + std::to_string(faceIndex) + " is out of range for face list (size " +
                std::to_string(allFaces.size()) + ").");
            }

            const Face& face = allFaces[faceIndex];

            // Ensure the face properties were calculated
            if(!face.geometricPropertiesCalculated) {
                throw std::runtime_error("Error in Cell " + std::to_string(id) + " calculation: " +
                                         "Geometric properties for bounding Face " + std::to_string(face.id) + " were not calculated.");
            }

            // Calculate cell centroid and volume
            Vector adjustedNormal = faceSigns[i] * face.normal;

            Scalar cf_dot_Sf = faceSigns[i] * face.area * dot(face.centroid, face.normal);

            centroidSum.x += adjustedNormal.x * face.x2_integral;
            centroidSum.y += adjustedNormal.y * face.y2_integral;
            centroidSum.z += adjustedNormal.z * face.z2_integral;

            volume += cf_dot_Sf;
        }

        volume /= S(3.0);
        if (std::abs(volume) > DIVISION_TOLERANCE) {
            centroid = centroidSum / (S(2.0) * volume);
        } else {
            throw std::runtime_error("Cell " + std::to_string(id) + " has zero volume");
        }

        if (volume < S(0.0)) {
            throw std::runtime_error("Error: Cell " + std::to_string(id) + " calculated negative volume (" + std::to_string(volume) + "). Check face normal conventions and mesh connectivity.");
        }
        else if (std::abs(volume) < VOLUME_TOLERANCE) {
            throw std::runtime_error("Warning: Cell " + std::to_string(id) + " calculated volume is near-zero (" + std::to_string(volume) + "). Check mesh quality or face properties.");
        }
        else {
            geometricPropertiesCalculated = true;
        }
    }
};

/* This function is used to print the cell to the console
 * It prints the cell's id, faces, neighbours, and calculated properties if available
 * The function sets the precision for floating Vector output within this scope
 * The function then prints the cell's id, faces, neighbours, and calculated properties if available
 * The function then restores the precision for floating Vector output within this scope
 * The function then returns the ostream object
 */
inline std::ostream& operator<<(std::ostream& os, const Cell& c) {
    os << "Cell(ID: " << c.id
       << ", Faces: [";
    // Print face indices
    for (size_t i = 0; i < c.faceIndices.size(); ++i) {
        os << c.faceIndices[i] << (i == c.faceIndices.size() - 1 ? "" : ", ");
    }
    os << "], Neighbours: [";
    // Print neighbour indices
     for (size_t i = 0; i < c.neighbourCellIndices.size(); ++i) {
        os << c.neighbourCellIndices[i] << (i == c.neighbourCellIndices.size() - 1 ? "" : ", ");
    }
    os << "]";

    // Print calculated properties if available
    if (c.geometricPropertiesCalculated) {
        std::ios_base::fmtflags flags = os.flags(); 
        int                      prec = os.precision(); 
        os << std::fixed << std::setprecision(6); 
        os << ", Volume: " << c.volume
           << ", Centroid: " << c.centroid;
        os.flags(flags);
        os.precision(prec);
    } else {
        os << ", Geometry: N/A";
    }
    os << ")";
    return os;
}

#endif
