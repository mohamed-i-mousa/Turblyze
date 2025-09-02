#include "Cell.h"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>

void Cell::calculateGeometricProperties(const std::vector<Face>& allFaces)
{
    geometricPropertiesCalculated = false;
    volume = 0.0;
    centroid = Vector(0.0, 0.0, 0.0);
    Vector centroidSum(0.0, 0.0, 0.0);

    for (size_t i = 0; i < faceIndices.size(); ++i)
    { 
        size_t faceIndex = faceIndices[i];
        const Face& face = allFaces[faceIndex];

        if(!face.geometricPropertiesCalculated) 
        {
            throw std::runtime_error
                (
                    "Error in Cell " + std::to_string(id) + " calculation:"
                  + " Geometric properties for bounding Face "
                  + std::to_string(face.id) + " were not calculated."
                );
        }
        
        Vector adjustedNormal = faceSigns[i] * face.normal;

        Scalar cf_dot_Sf = 
            faceSigns[i] * face.area * dot(face.centroid, face.normal);

        centroidSum.x += adjustedNormal.x * face.x2_integral;
        centroidSum.y += adjustedNormal.y * face.y2_integral;
        centroidSum.z += adjustedNormal.z * face.z2_integral;

        volume += cf_dot_Sf;
    }

    volume /= S(3.0);
    
    if (std::abs(volume) > smallValue) 
    {
        centroid = centroidSum / (S(2.0) * volume);
    } 
    else 
    {
        throw std::runtime_error
        (
            "Cell " + std::to_string(id) + " has zero volume"
        );
    }

    if (volume < S(0.0)) 
    {
        throw std::runtime_error
        (
            "Error: Cell " + std::to_string(id) 
          + " calculated negative volume (" + std::to_string(volume)
          + "). Check face normal conventions and mesh connectivity."
        );
    }
    else
    {
        geometricPropertiesCalculated = true;
    }
}

std::ostream& operator<<(std::ostream& os, const Cell& c)
{
    os  << "Cell(ID: " << c.id << ", Faces: [";
    
    for (size_t i = 0; i < c.faceIndices.size(); ++i) 
    {
        os  << c.faceIndices[i] << (i == c.faceIndices.size() - 1 ? "" : ", ");
    }

    os  << "], Neighbours: [";
    
    for (size_t i = 0; i < c.neighbourCellIndices.size(); ++i)
    {
        os << c.neighbourCellIndices[i] << (i == c.neighbourCellIndices.size() - 1 ? "" : ", ");
    }
    
    os  << "]";

    if (c.geometricPropertiesCalculated)
    {
        std::ios_base::fmtflags flags = os.flags(); 
        int prec = os.precision(); 

        os  << std::fixed << std::setprecision(6); 
        os  << ", Volume: " << c.volume
            << ", Centroid: " << c.centroid;
            
        os.flags(flags);
        os.precision(prec);
    }
    else
    {
        os  << ", Geometry: N/A";
    }
    
    os << ")";
    
    return os;
}
