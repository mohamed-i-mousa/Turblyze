/******************************************************************************
 * @file Face.cpp
 * @brief Implementation of face geometric properties and operations
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include <ostream>
#include <iomanip>

#include "Face.hpp"
#include "Cell.hpp"

// ********************** Geometric Property Methods **********************

void Face::calculateGeometricProperties(std::span<const Vector> allNodes)
{
    geometricPropertiesCalculated_ = false;
    const size_t numNodes = nodeIndices_.size();

    for (size_t nodeIdx : nodeIndices_)
    {
        if (nodeIdx >= allNodes.size())
        {
            throw
                std::out_of_range
                (
                    "Error calculating properties for Face "
                  + std::to_string(idx_) + ": Node index "
                  + std::to_string(nodeIdx)
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
        Scalar x2Formula =
            p1.x()*p1.x() + p2.x()*p2.x() + p3.x()*p3.x()
          + p1.x()*p2.x() + p1.x()*p3.x() + p2.x()*p3.x();

        Scalar y2Formula =
            p1.y()*p1.y() + p2.y()*p2.y() + p3.y()*p3.y()
          + p1.y()*p2.y() + p1.y()*p3.y() + p2.y()*p3.y();

        Scalar z2Formula =
            p1.z()*p1.z() + p2.z()*p2.z() + p3.z()*p3.z()
          + p1.z()*p2.z() + p1.z()*p3.z() + p2.z()*p3.z();

        x2Integral_ = crossProd.x() * x2Formula / S(12.0);
        y2Integral_ = crossProd.y() * y2Formula / S(12.0);
        z2Integral_ = crossProd.z() * z2Formula / S(12.0);

        volumeContribution_ = dot(centroid_, crossProd) / S(2.0);

        geometricPropertiesCalculated_ = true;
    }
    // CASE 2: Face is "Polygon" (numNodes > 3)
    else
    {
        Vector faceCenter(0.0, 0.0, 0.0);

        for (size_t nodeIdx : nodeIndices_)
        {
            faceCenter += allNodes[nodeIdx];
        }

        faceCenter /= S(numNodes);

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
            const Vector& pCurr = allNodes[nodeIndices_[i]];
            const Vector& pNext = allNodes[nodeIndices_[(i + 1) % numNodes]];

            const Vector& p1Tri = faceCenter;
            const Vector& p2Tri = pCurr;
            const Vector& p3Tri = pNext;

            Vector vecATri = p2Tri - p1Tri;
            Vector vecBTri = p3Tri - p1Tri;

            Vector crossProdTri = cross(vecBTri, vecATri);
            Scalar triangleArea = S(0.5) * crossProdTri.magnitude();
            normalSum += crossProdTri;

            // Weighted integrals using each sub-triangle's normal
            // For divergence theorem: contribution = N_x * x2Formula / 12
            Scalar x2Formula =
                p1Tri.x() * p1Tri.x() + p2Tri.x() * p2Tri.x()
              + p3Tri.x() * p3Tri.x() + p1Tri.x() * p2Tri.x()
              + p1Tri.x() * p3Tri.x() + p2Tri.x() * p3Tri.x();

            Scalar y2Formula =
                p1Tri.y() * p1Tri.y() + p2Tri.y() * p2Tri.y()
              + p3Tri.y() * p3Tri.y() + p1Tri.y() * p2Tri.y()
              + p1Tri.y() * p3Tri.y() + p2Tri.y() * p3Tri.y();

            Scalar z2Formula =
                p1Tri.z() * p1Tri.z() + p2Tri.z() * p2Tri.z()
              + p3Tri.z() * p3Tri.z() + p1Tri.z() * p2Tri.z()
              + p1Tri.z() * p3Tri.z() + p2Tri.z() * p3Tri.z();

            x2Integral_ += crossProdTri.x() * x2Formula / S(12.0);
            y2Integral_ += crossProdTri.y() * y2Formula / S(12.0);
            z2Integral_ += crossProdTri.z() * z2Formula / S(12.0);

            Vector triangleCentroid = (p1Tri + p2Tri + p3Tri) / S(3.0);

            // Volume contribution for this sub-triangle:
            volumeContribution_ +=
                dot(triangleCentroid, crossProdTri) / S(2.0);

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
    os  << "Face(ID: " << f.idx() << ", Nodes: [";

    const auto& nodes = f.nodeIndices();
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        os  << nodes[i]
            << (i == nodes.size() - 1 ? "" : ", ");
    }

    os  <<  "], Owner: " << f.ownerCell() << ", Neighbor: "
        <<  (
                f.isBoundary() ? "Boundary"
              : std::to_string(f.neighborCell().value_or(0))
            );

    if (f.geometricPropertiesCalculated())
    {
        std::ios_base::fmtflags flags = os.flags();
        auto prec = os.precision();

        os  << std::fixed << std::setprecision(6);

        os  << ", Centroid: " << f.centroid()
            << ", Area: "   << f.projectedArea()
            << ", Normal: " << f.normal();

        os  << ", dPfMag: " << f.dPfMag();

        if (f.dNfMag().has_value())
        {
            os  << ", dNfMag: " << f.dNfMag().value();
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
