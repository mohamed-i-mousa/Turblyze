#ifndef CELL_H
#define CELL_H

#include <vector>
#include <stdexcept>  // for runtime_error, out_of_range
#include <cmath>
#include <iostream>   // for std::cerr, std::ostream
#include <iomanip>    // for std::setprecision, std::fixed in operator<<

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"


struct Cell {
    // ----- Members ----- //
    size_t id = 0;
    std::vector<size_t> faceIndices;
    std::vector<size_t> neighbourCellIndices;
    std::vector<int> faceSigns;                     // This to adjust normal vector direction to be always pointing outward 
                                                    // If the cell is the owner (sign = 1), if the cell is neighbour (sign = -1)
    Vector centroid;
    Scalar volume = 0.0;
    bool geometricPropertiesCalculated = false;

    // ----- Constructors ----- //

    // Default constructor 
    Cell() = default;

    // Parameterized constructor 
    Cell(size_t cellId, const std::vector<size_t>& faces, const std::vector<size_t>& neighbours, const std::vector<int>& signs)
        : id(cellId),
          faceIndices(faces),
          faceSigns(signs),
          neighbourCellIndices(neighbours)
          {}

    // ----- Member Methods ----- //
    
    void calculateGeometricProperties(const std::vector<Face>& allFaces) {
        geometricPropertiesCalculated = false;
        volume = 0.0;
        centroid = Vector(0.0, 0.0, 0.0);
        Vector centroidSum(0.0, 0.0, 0.0);

        if (faceIndices.empty()) {
            std::cerr << "Warning: Cell " << id << "has no assigned faces. Calculations cannot proceed!";
            return;
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
            Vector Sf;
            Vector adjustedNormal = faceSigns[i] * face.normal;

            Sf = adjustedNormal * face.area;

            Scalar cf_dot_Sf = dot(face.centroid, Sf);

            centroidSum.x += adjustedNormal.x * face.x2_integral;
            centroidSum.y += adjustedNormal.y * face.y2_integral;
            centroidSum.z += adjustedNormal.z * face.z2_integral;

            volume += cf_dot_Sf;
        }

        volume /= 3.0;
        centroid = centroidSum / (2.0 * volume);

        // Check for invalid volume
        if (std::fabs(volume) < 1e-12) {
            throw std::runtime_error("Warning: Cell " + std::to_string(id) + " calculated volume is near-zero (" + std::to_string(volume) + "). Check mesh quality or face properties. Centroid calculation skipped.");

        }
        // Check for negative volumes
        else if (volume < 0.0) {
            throw std::runtime_error("Error: Cell " + std::to_string(id) + " calculated negative volume (" + std::to_string(volume) + "). Check face normal conventions and mesh connectivity. Centroid calculation skipped.");
        }
        else {
            geometricPropertiesCalculated = true;
        }
    }
};

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
        std::ios_base::fmtflags flags = os.flags(); // Save format flags
        int                      prec = os.precision(); // Save precision
        os << std::fixed << std::setprecision(6); // Apply formatting
        os << ", Volume: " << c.volume
           << ", Centroid: " << c.centroid;
        os.flags(flags); // Restore flags
        os.precision(prec); // Restore precision
    } else {
        os << ", Geometry: N/A"; // Indicate properties not calculated or failed
    }
    os << ")"; // Closing parenthesis
    return os;
}

#endif
