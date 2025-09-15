/******************************************************************************
 * @file Face.cpp
 * @brief Implementation of face geometric properties and operations
 *****************************************************************************/

#include "Face.hpp"
#include "Cell.hpp"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>

void Face::calculateGeometricProperties(const std::vector<Vector>& allNodes) 
{
    geometricPropertiesCalculated_ = false;         
    const size_t nNodes = nodeIndices_.size();

    for (size_t i = 0; i < nNodes; ++i) 
    {
        if (nodeIndices_[i] >= allNodes.size()) 
        {
            throw std::out_of_range
                (
                    "Error calculating properties for Face "
                  + std::to_string(id_) + ": Node index "
                  + std::to_string(nodeIndices_[i]) + " out of range ("
                  + "Node list size: " + std::to_string(allNodes.size())
                  + ")."
                );
        }
    }

    // CASE 1: Face is "Triangle" (nNodes == 3)
    if (nNodes == 3) 
    {
        const Vector& p1 = allNodes[nodeIndices_[0]];
        const Vector& p2 = allNodes[nodeIndices_[1]];
        const Vector& p3 = allNodes[nodeIndices_[2]];

        centroid_ = (p1 + p2 + p3) / S(3.0);

        Vector vecA = p2 - p1;
        Vector vecB = p3 - p1;

        Vector crossProd = cross(vecB, vecA);
        Scalar crossProdMag = crossProd.magnitude();

        area_ = S(0.5) * crossProdMag;
        normal_ = crossProd / (crossProdMag + vSmallValue);

        x2_integral_ = 
            (p1.x()*p1.x() + p2.x()*p2.x() + p3.x()*p3.x() 
            + p1.x()*p2.x() + p1.x()*p3.x() + p2.x()*p3.x())
            * area_ / S(6.0);
        y2_integral_ = 
            (p1.y()*p1.y() + p2.y()*p2.y() + p3.y()*p3.y()
            + p1.y()*p2.y() + p1.y()*p3.y() + p2.y()*p3.y())
            * area_ / S(6.0);
        z2_integral_ = 
            (p1.z()*p1.z() + p2.z()*p2.z() + p3.z()*p3.z()
            + p1.z()*p2.z() + p1.z()*p3.z() + p2.z()*p3.z())
            * area_ / S(6.0);

        geometricPropertiesCalculated_ = true;
    }
    // CASE 2: Face is "Polygon" (nNodes > 3)
    else 
    {
        Vector faceCenter(0.0, 0.0, 0.0);

        for (size_t i = 0; i < nNodes; ++i) 
        {
            faceCenter += allNodes[nodeIndices_[i]];
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
            const Vector& p_i = allNodes[nodeIndices_[i]];
            // Handles the end point
            const Vector& p_next = allNodes[nodeIndices_[(i + 1) % nNodes]];
            
            const Vector& p1_tri = faceCenter;
            const Vector& p2_tri = p_i;
            const Vector& p3_tri = p_next;

            Vector vecA_tri = p2_tri - p1_tri;
            Vector vecB_tri = p3_tri - p1_tri;

            Vector crossProd_tri = cross(vecB_tri, vecA_tri);
            Scalar triangleArea = S(0.5) * crossProd_tri.magnitude();
            normalSum += crossProd_tri;

            Scalar x2_part = 
                (p1_tri.x() * p1_tri.x() + p2_tri.x() * p2_tri.x()
               + p3_tri.x() * p3_tri.x() + p1_tri.x() * p2_tri.x()
               + p1_tri.x() * p3_tri.x() + p2_tri.x() * p3_tri.x())
               * triangleArea / S(6.0);
            Scalar y2_part = 
                (p1_tri.y() * p1_tri.y() + p2_tri.y() * p2_tri.y()
               + p3_tri.y() * p3_tri.y() + p1_tri.y() * p2_tri.y()
               + p1_tri.y() * p3_tri.y() + p2_tri.y() * p3_tri.y())
               * triangleArea / S(6.0);
            Scalar z2_part = 
                (p1_tri.z() * p1_tri.z() + p2_tri.z() * p2_tri.z()
               + p3_tri.z() * p3_tri.z() + p1_tri.z() * p2_tri.z()
               + p1_tri.z() * p3_tri.z() + p2_tri.z() * p3_tri.z())
               * triangleArea / S(6.0);

            x2_integral_ += x2_part;
            y2_integral_ += y2_part;
            z2_integral_ += z2_part;

            Vector triangleCentroid = (p1_tri + p2_tri + p3_tri) / S(3.0);
            totalArea += triangleArea;
            weightedCentroidSum += triangleCentroid * triangleArea;
        }

        area_ = totalArea;

        centroid_ = weightedCentroidSum / (area_ + vSmallValue);

        normal_ = normalSum.normalized();

        geometricPropertiesCalculated_ = true;
    }
}

std::ostream& operator<<(std::ostream& os, const Face& f)
{
   os  << "Face(ID: " << f.id_
       << ", Nodes: [";
   
   for (size_t i = 0; i < f.nodeIndices_.size(); ++i)
   {
       os  << f.nodeIndices_[i] << (i == f.nodeIndices_.size() - 1 ? "" : ", ");
   }
   os  << "], Owner: " << f.ownerCell_
       << ", Neighbor: " << (
                               f.isBoundary() ? "Boundary"
                               : std::to_string(f.neighborCell_.value_or(0))
                             );

   if (f.geometricPropertiesCalculated_) 
   {
       std::ios_base::fmtflags flags = os.flags();
       int prec = os.precision();

       os  << std::fixed << std::setprecision(6);
       os  << ", Centroid: " << f.centroid_
           << ", Area: " << f.area_
           << ", Normal: " << f.normal_;
       
       if (f.distancePropertiesCalculated_) 
       {
           os  << ", d_Pf_mag: " << f.d_Pf_mag_;

           if (f.d_Nf_mag_.has_value()) 
           {
               os << ", d_Nf_mag: " << f.d_Nf_mag_.value();
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
    d_Pf_ = centroid_ - allCells[ownerCell_].centroid();
    d_Pf_mag_ = d_Pf_.magnitude();
    
    e_Pf_ = d_Pf_ / (d_Pf_mag_ + vSmallValue);

    // Calculate d_Nf only for internal faces
    if (!isBoundary())
    {
        size_t N = neighborCell_.value();
        
        Vector d_Nf_vec = centroid_ - allCells[N].centroid();
        Scalar d_Nf_magnitude = d_Nf_vec.magnitude();

        d_Nf_ = d_Nf_vec;
        d_Nf_mag_ = d_Nf_magnitude;
        e_Nf_ = d_Nf_vec / (d_Nf_magnitude + vSmallValue);
    }
    distancePropertiesCalculated_ = true;
}

// Explicit template instantiation for the types used in the codebase
template void Face::calculateDistanceProperties(const std::vector<Cell>& allCells);
