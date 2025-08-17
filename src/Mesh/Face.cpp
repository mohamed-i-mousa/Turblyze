#include "Face.h"
#include "Cell.h"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>

void Face::calculateGeometricProperties(const std::vector<Vector>& allNodes) 
{
    geometricPropertiesCalculated = false;         
    const size_t nNodes = nodeIndices.size();

    for (size_t i = 0; i < nNodes; ++i) 
    {
        if (nodeIndices[i] >= allNodes.size()) 
        {
            throw std::out_of_range
                (
                    "Error calculating properties for Face "
                  + std::to_string(id) + ": Node index "
                  + std::to_string(nodeIndices[i]) + " out of range ("
                  + "Node list size: " + std::to_string(allNodes.size())
                  + ")."
                );
        }
    }

    // CASE 1: Face is "Triangle" (nNodes == 3)
    if (nNodes == 3) 
    {
        const Vector& p1 = allNodes[nodeIndices[0]];
        const Vector& p2 = allNodes[nodeIndices[1]];
        const Vector& p3 = allNodes[nodeIndices[2]];

        centroid = (p1 + p2 + p3) / S(3.0);

        Vector vecA = p2 - p1;
        Vector vecB = p3 - p1;

        Vector crossProd = cross(vecB, vecA);
        Scalar crossProdMag = crossProd.magnitude();

        if (std::abs(crossProdMag) < AREA_TOLERANCE) 
        {
            throw std::runtime_error
                (
                    "Face " + std::to_string(id) + " is degenerate."
                );
        } 
        else 
        {
            area = S(0.5) * crossProdMag;
            normal = crossProd / crossProdMag;

            x2_integral = 
                (p1.x*p1.x + p2.x*p2.x + p3.x*p3.x 
               + p1.x*p2.x + p1.x*p3.x + p2.x*p3.x)
               * area / S(6.0);
            y2_integral = 
                (p1.y*p1.y + p2.y*p2.y + p3.y*p3.y
               + p1.y*p2.y + p1.y*p3.y + p2.y*p3.y)
               * area / S(6.0);
            z2_integral = 
                (p1.z*p1.z + p2.z*p2.z + p3.z*p3.z
               + p1.z*p2.z + p1.z*p3.z + p2.z*p3.z)
               * area / S(6.0);

            geometricPropertiesCalculated = true;
        }
    }
    // CASE 2: Face is "Polygon" (nNodes > 3)
    else 
    {
        Vector faceCenter(0.0, 0.0, 0.0);

        for (size_t i = 0; i < nNodes; ++i) 
        {
            faceCenter += allNodes[nodeIndices[i]];
        }

        if (nNodes > 0) 
        {
            faceCenter /= S(nNodes);
        }

        Scalar totalArea = 0.0;
        Vector weightedCentroidSum(0.0, 0.0, 0.0);
        Vector normalSum(0.0, 0.0, 0.0);

        for (size_t i = 0; i < nNodes; ++i) 
        {
            const Vector& p_i = allNodes[nodeIndices[i]];
            // Handles the end point
            const Vector& p_next = allNodes[nodeIndices[(i + 1) % nNodes]];
            
            const Vector& p1_tri = faceCenter;
            const Vector& p2_tri = p_i;
            const Vector& p3_tri = p_next;

            Vector vecA_tri = p2_tri - p1_tri;
            Vector vecB_tri = p3_tri - p1_tri;

            Vector crossProd_tri = cross(vecB_tri, vecA_tri);
            Scalar triangleArea = S(0.5) * crossProd_tri.magnitude();
            normalSum += crossProd_tri;

            Scalar x2_part = 
                (p1_tri.x * p1_tri.x + p2_tri.x * p2_tri.x
               + p3_tri.x * p3_tri.x + p1_tri.x * p2_tri.x
               + p1_tri.x * p3_tri.x + p2_tri.x * p3_tri.x)
               * triangleArea / S(6.0);
            Scalar y2_part = 
                (p1_tri.y * p1_tri.y + p2_tri.y * p2_tri.y
               + p3_tri.y * p3_tri.y + p1_tri.y * p2_tri.y
               + p1_tri.y * p3_tri.y + p2_tri.y * p3_tri.y)
               * triangleArea / S(6.0);
            Scalar z2_part = 
                (p1_tri.z * p1_tri.z + p2_tri.z * p2_tri.z
               + p3_tri.z * p3_tri.z + p1_tri.z * p2_tri.z
               + p1_tri.z * p3_tri.z + p2_tri.z * p3_tri.z)
               * triangleArea / S(6.0);

            x2_integral += x2_part;
            y2_integral += y2_part;
            z2_integral += z2_part;

            if (triangleArea > AREA_TOLERANCE)
            {
                Vector triangleCentroid = (p1_tri + p2_tri + p3_tri) / S(3.0);
                totalArea += triangleArea;
                weightedCentroidSum += triangleCentroid * triangleArea;
            } 
            else 
            {
                throw std::runtime_error
                    (
                        "Warning: Polygonal Face " + std::to_string(id)
                      + " has near-zero total area."
                    );
            }
        }

        area = totalArea;

        if (area < AREA_TOLERANCE)
        {
             throw std::runtime_error
                (
                    "Warning: Polygonal Face " + std::to_string(id)
                  + " has near-zero total area. Setting area=0, "
                  + "centroid/normal=(0,0,0)."
                );
        } 
        else 
        {
             if (std::abs(area) > DIVISION_TOLERANCE)
             {
                 centroid = weightedCentroidSum / area;
             } 
             else 
             {
                 throw std::runtime_error
                    (
                        "Warning: Polygonal Face " + std::to_string(id)
                      + " has near-zero total area."
                    );
             }
             normal = normalSum.normalized();
             geometricPropertiesCalculated = true;
        }
    }
}

std::ostream& operator<<(std::ostream& os, const Face& f)
{
   os  << "Face(ID: " << f.id
       << ", Nodes: [";
   
   for (size_t i = 0; i < f.nodeIndices.size(); ++i)
   {
       os  << f.nodeIndices[i] << (i == f.nodeIndices.size() - 1 ? "" : ", ");
   }
   os  << "], Owner: " << f.ownerCell
       << ", Neighbour: " << (
                               f.isBoundary() ? "Boundary"
                               : std::to_string(f.neighbourCell.value_or(0))
                             );

   if (f.geometricPropertiesCalculated) 
   {
       std::ios_base::fmtflags flags = os.flags();
       int prec = os.precision();

       os  << std::fixed << std::setprecision(6);
       os  << ", Centroid: " << f.centroid
           << ", Area: " << f.area
           << ", Normal: " << f.normal;
       
       if (f.distancePropertiesCalculated) 
       {
           os  << ", d_Pf_mag: " << f.d_Pf_mag;

           if (f.d_Nf_mag.has_value()) 
           {
               os << ", d_Nf_mag: " << f.d_Nf_mag.value();
           }
       }
       
       os.flags(flags);
       os.precision(prec);
   }
   else
   {
       os  << ", Geometry: N/A";
   }
   os  << ")";

   return os;
}

template<typename CellContainer>
void Face::calculateDistanceProperties(const CellContainer& allCells)
{ 
    d_Pf = centroid - allCells[ownerCell].centroid;
    d_Pf_mag = d_Pf.magnitude();
    
    if (d_Pf_mag > DIVISION_TOLERANCE)
    {
        e_Pf = d_Pf / d_Pf_mag;
    } 
    else 
    {
        throw std::runtime_error
            (
                "Face " + std::to_string(id)
              + ": distance from owner cell to face is nearly zero."
            );
    }

    // Calculate d_Nf only for internal faces
    if (!isBoundary())
    {
        size_t N = neighbourCell.value();
        
        Vector d_Nf_vec = centroid - allCells[N].centroid;
        Scalar d_Nf_magnitude = d_Nf_vec.magnitude();
        
        if (d_Nf_magnitude > DIVISION_TOLERANCE)
        {
            d_Nf = d_Nf_vec;
            d_Nf_mag = d_Nf_magnitude;
            e_Nf = d_Nf_vec / d_Nf_magnitude;
        } 
        else 
        {
            throw std::runtime_error
                (
                    "Face " + std::to_string(id)
                  + ": distance from neighbor cell to face is nearly zero."
                );
        }
    }
    distancePropertiesCalculated = true;
}

// Explicit template instantiation for the types used in the codebase
template void Face::calculateDistanceProperties(const std::vector<Cell>& allCells);
