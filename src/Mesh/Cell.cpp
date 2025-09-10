/******************************************************************************
 * @file Cell.cpp
 * @brief Implementation of cell geometric properties and operations
 *****************************************************************************/

#include "Cell.h"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>

void Cell::calculateGeometricProperties(const std::vector<Face>& allFaces)
{
    geometricPropertiesCalculated_ = false;
    volume_ = 0.0;
    centroid_ = Vector(0.0, 0.0, 0.0);
    Vector centroidSum(0.0, 0.0, 0.0);

    for (size_t i = 0; i < faceIndices_.size(); ++i)
    { 
        size_t faceIndex = faceIndices_[i];
        const Face& face = allFaces[faceIndex];

        if(!face.geometricPropertiesCalculated()) 
        {
            throw std::runtime_error
                (
                    "Error in Cell " + std::to_string(id_) + " calculation:"
                  + " Geometric properties for bounding Face "
                  + std::to_string(face.id()) + " were not calculated."
                );
        }
        
        Vector adjustedNormal = faceSigns_[i] * face.normal();

        Scalar cf_dot_Sf = 
            faceSigns_[i] * face.area() * dot(face.centroid(), face.normal());

        centroidSum.setX(centroidSum.x() + adjustedNormal.x() * face.x2_integral());
        centroidSum.setY(centroidSum.y() + adjustedNormal.y() * face.y2_integral());
        centroidSum.setZ(centroidSum.z() + adjustedNormal.z() * face.z2_integral());

        volume_ += cf_dot_Sf;
    }

    volume_ /= S(3.0);
    
    if (std::abs(volume_) > smallValue) 
    {
        centroid_ = centroidSum / (S(2.0) * volume_);
    } 
    else 
    {
        throw std::runtime_error
        (
            "Cell " + std::to_string(id_) + " has zero volume"
        );
    }

    if (volume_ < S(0.0)) 
    {
        throw std::runtime_error
        (
            "Error: Cell " + std::to_string(id_) 
          + " calculated negative volume (" + std::to_string(volume_)
          + "). Check face normal conventions and mesh connectivity."
        );
    }
    else
    {
        geometricPropertiesCalculated_ = true;
    }
}

std::ostream& operator<<(std::ostream& os, const Cell& c)
{
    os  << "Cell(ID: " << c.id_ << ", Faces: [";
    
    for (size_t i = 0; i < c.faceIndices_.size(); ++i) 
    {
        os  << c.faceIndices_[i] << (i == c.faceIndices_.size() - 1 ? "" : ", ");
    }

    os  << "], Neighbors: [";
    
    for (size_t i = 0; i < c.neighborCellIndices_.size(); ++i)
    {
        os << c.neighborCellIndices_[i] << (i == c.neighborCellIndices_.size() - 1 ? "" : ", ");
    }
    
    os  << "]";

    if (c.geometricPropertiesCalculated_)
    {
        std::ios_base::fmtflags flags = os.flags(); 
        int prec = os.precision(); 

        os  << std::fixed << std::setprecision(6); 
        os  << ", Volume: " << c.volume_
            << ", Centroid: " << c.centroid_;
            
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
