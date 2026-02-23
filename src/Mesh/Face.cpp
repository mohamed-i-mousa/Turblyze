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

// ********************** Geometric Property Methods **********************

void Face::calculateGeometricProperties(const std::vector<Vector>& allNodes)
{
    geometricPropertiesCalculated_ = false;
    const size_t numNodes = nodeIndices_.size();

    for (size_t i = 0; i < numNodes; ++i)
    {
        if (nodeIndices_[i] >= allNodes.size())
        {
            throw
                std::out_of_range
                (
                    "Error calculating properties for Face "
                  + std::to_string(idx_) + ": Node index "
                  + std::to_string(nodeIndices_[i])
                  + " out of range ("
                  + "Node list size: "
                  + std::to_string(allNodes.size())
                  + ")."
                );
        }
    }

    // CASE 1: Face is "Triangle" (numNodes == 3)
    if (numNodes == 3)
    {
        const Vector& p1 = allNodes[nodeIndices_[0]];
        const Vector& p2 = allNodes[nodeIndices_[1]];
        const Vector& p3 = allNodes[nodeIndices_[2]];

        centroid_ = (p1 + p2 + p3) / S(3.0);

        Vector vecA = p2 - p1;
        Vector vecB = p3 - p1;

        Vector crossProd = cross(vecB, vecA);
        Scalar crossProdMag = crossProd.magnitude();

        projectedArea_ = S(0.5) * crossProdMag;

        // For planar triangles, contact = projected
        contactArea_ = projectedArea_;

        normal_ = crossProd / (crossProdMag + vSmallValue);

        // Second moment integrals weighted by normal component
        Scalar x2_formula =
            p1.x()*p1.x() + p2.x()*p2.x() + p3.x()*p3.x()
          + p1.x()*p2.x() + p1.x()*p3.x() + p2.x()*p3.x();

        Scalar y2_formula =
            p1.y()*p1.y() + p2.y()*p2.y() + p3.y()*p3.y()
          + p1.y()*p2.y() + p1.y()*p3.y() + p2.y()*p3.y();

        Scalar z2_formula =
            p1.z()*p1.z() + p2.z()*p2.z() + p3.z()*p3.z()
          + p1.z()*p2.z() + p1.z()*p3.z() + p2.z()*p3.z();

        x2Integral_ = crossProd.x() * x2_formula / S(12.0);
        y2Integral_ = crossProd.y() * y2_formula / S(12.0);
        z2Integral_ = crossProd.z() * z2_formula / S(12.0);

        volumeContribution_ = dot(centroid_, crossProd) / S(2.0);

        geometricPropertiesCalculated_ = true;
    }
    // CASE 2: Face is "Polygon" (numNodes > 3)
    else
    {
        Vector faceCenter(0.0, 0.0, 0.0);

        for (size_t i = 0; i < numNodes; ++i)
        {
            faceCenter += allNodes[nodeIndices_[i]];
        }

        if (numNodes > 0)
        {
            faceCenter /= S(numNodes);
        }

        Vector weightedCentroidSum(0.0, 0.0, 0.0);
        Vector normalSum(0.0, 0.0, 0.0);
        Scalar weightedAreaSum = 0.0;

        // Reset integrals (for accumulation)
        x2Integral_ = 0.0;
        y2Integral_ = 0.0;
        z2Integral_ = 0.0;
        volumeContribution_ = 0.0;

        for (size_t i = 0; i < numNodes; ++i)
        {
            const Vector& p_i = allNodes[nodeIndices_[i]];
            const Vector& p_next = allNodes[nodeIndices_[(i + 1) % numNodes]];

            const Vector& p1_tri = faceCenter;
            const Vector& p2_tri = p_i;
            const Vector& p3_tri = p_next;

            Vector vecA_tri = p2_tri - p1_tri;
            Vector vecB_tri = p3_tri - p1_tri;

            Vector crossProd_tri = cross(vecB_tri, vecA_tri);
            Scalar triangleArea = S(0.5) * crossProd_tri.magnitude();
            normalSum += crossProd_tri;

            // Weighted integrals using each sub-triangle's normal
            // For divergence theorem: contribution = N_x * x2_formula / 12
            Scalar x2_formula =
                p1_tri.x() * p1_tri.x() + p2_tri.x() * p2_tri.x()
              + p3_tri.x() * p3_tri.x() + p1_tri.x() * p2_tri.x()
              + p1_tri.x() * p3_tri.x() + p2_tri.x() * p3_tri.x();

            Scalar y2_formula =
                p1_tri.y() * p1_tri.y() + p2_tri.y() * p2_tri.y()
              + p3_tri.y() * p3_tri.y() + p1_tri.y() * p2_tri.y()
              + p1_tri.y() * p3_tri.y() + p2_tri.y() * p3_tri.y();

            Scalar z2_formula =
                p1_tri.z() * p1_tri.z() + p2_tri.z() * p2_tri.z()
              + p3_tri.z() * p3_tri.z() + p1_tri.z() * p2_tri.z()
              + p1_tri.z() * p3_tri.z() + p2_tri.z() * p3_tri.z();

            x2Integral_ += crossProd_tri.x() * x2_formula / S(12.0);
            y2Integral_ += crossProd_tri.y() * y2_formula / S(12.0);
            z2Integral_ += crossProd_tri.z() * z2_formula / S(12.0);

            Vector triangleCentroid = (p1_tri + p2_tri + p3_tri) / S(3.0);

            // Volume contribution for this sub-triangle:
            volumeContribution_ +=
                dot(triangleCentroid, crossProd_tri) / S(2.0);

            weightedAreaSum += triangleArea;
            weightedCentroidSum += triangleCentroid * triangleArea;
        }

        projectedArea_ = normalSum.magnitude() / S(2.0);

        contactArea_ = weightedAreaSum;

        // Centroid uses contact area weighting (sum of triangle areas)
        centroid_ = weightedCentroidSum / (weightedAreaSum + vSmallValue);

        normal_ = normalSum.normalized();

        geometricPropertiesCalculated_ = true;
    }
}

std::ostream& operator<<(std::ostream& os, const Face& f)
{
    os  << "Face(ID: " << f.idx_ << ", Nodes: [";

    for (size_t i = 0; i < f.nodeIndices_.size(); ++i)
    {
        os  << f.nodeIndices_[i]
            << (i == f.nodeIndices_.size() - 1 ? "" : ", ");
    }

    os  <<  "], Owner: " << f.ownerCell_ << ", Neighbor: "
        <<  (
                f.isBoundary() ? "Boundary"
              : std::to_string(f.neighborCell_.value_or(0))
            );

    if (f.geometricPropertiesCalculated_)
    {
        std::ios_base::fmtflags flags = os.flags();
        int prec = os.precision();

        os  << std::fixed << std::setprecision(6);

        os  << ", Centroid: " << f.centroid_ 
            << ", Area: "   << f.projectedArea_
            << ", Normal: " << f.normal_;

        if (f.distancePropertiesCalculated_)
        {
            os  << ", dPfMag: " << f.dPfMag_;

            if (f.dNfMag_.has_value())
            {
                os  << ", dNfMag: " << f.dNfMag_.value();
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

// *********************** Distance Property Methods **********************

template<typename CellContainer>
void Face::calculateDistanceProperties(const CellContainer& allCells)
{ 
    dPf_ = centroid_ - allCells[ownerCell_].centroid();
    dPfMag_ = dPf_.magnitude();

    ePf_ = dPf_ / (dPfMag_ + vSmallValue);

    // Calculate dNf only for internal faces
    if (!isBoundary())
    {
        size_t N = neighborCell_.value();

        Vector dNfVec = centroid_ - allCells[N].centroid();
        Scalar dNfMagnitude = dNfVec.magnitude();

        dNf_ = dNfVec;
        dNfMag_ = dNfMagnitude;
        eNf_ = dNfVec / (dNfMagnitude + vSmallValue);
    }
    distancePropertiesCalculated_ = true;
}

// For header decoupling and to avoid circular dependencies
template void Face::calculateDistanceProperties
(
    const std::vector<Cell>& allCells
);