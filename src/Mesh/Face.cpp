/******************************************************************************
 * @file Face.cpp
 * @brief Implementation of face geometric properties and operations
 *****************************************************************************/

#include "Face.hpp"

#include <cmath>
#include <ostream>
#include <iomanip>

#include "ErrorHandler.hpp"


// ************************ Geometric Property Methods ************************

FaceIntegrals Face::calculateGeometricProperties
(
    std::span<const Vector> allNodes
)
{
    geometricPropertiesCalculated_ = false;
    const size_t numNodes = nodeIndices_.size();

    for (size_t nodeIdx : nodeIndices_)
    {
        if (nodeIdx >= allNodes.size())
        {
            FatalError
            (
                "Error calculating properties for Face "
              + std::to_string(idx_) + ": Node index "
              + std::to_string(nodeIdx)
              + " out of range (Node list size: "
              + std::to_string(allNodes.size()) + ")."
            );
        }
    }

    FaceIntegrals integrals;

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

        if (crossProdMag < vSmallValue)
        {
            FatalError
            (
                "Face " + std::to_string(idx_)
              + " is geometrically degenerate (zero area)."
            );
        }

        projectedArea_ = S(0.5) * crossProdMag;

        // For planar triangles, contact = projected
        contactArea_ = projectedArea_;

        normal_ = crossProd / crossProdMag;

        // Second moment integrals weighted by normal component
        integrals.x2 =
            crossProd.x() * secondMoment(p1.x(), p2.x(), p3.x()) / S(12.0);
        integrals.y2 =
            crossProd.y() * secondMoment(p1.y(), p2.y(), p3.y()) / S(12.0);
        integrals.z2 =
            crossProd.z() * secondMoment(p1.z(), p2.z(), p3.z()) / S(12.0);

        integrals.volume = dot(centroid_, crossProd) / S(2.0);

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

        Vector weightedCentroidSum{};
        Vector normalSum{};
        Scalar weightedAreaSum = 0.0;

        for (size_t nodeIdx = 0; nodeIdx < numNodes; ++nodeIdx)
        {
            const Vector& pCurr = allNodes[nodeIndices_[nodeIdx]];
            const Vector& pNext = allNodes[nodeIndices_[(nodeIdx + 1) % numNodes]];

            const Vector& p1Tri = faceCenter;
            const Vector& p2Tri = pCurr;
            const Vector& p3Tri = pNext;

            Vector vecATri = p2Tri - p1Tri;
            Vector vecBTri = p3Tri - p1Tri;

            Vector crossProdTri = cross(vecBTri, vecATri);
            Scalar triangleArea = S(0.5) * crossProdTri.magnitude();
            normalSum += crossProdTri;

            // Weighted integrals using each sub-triangle's normal
            integrals.x2 +=
                crossProdTri.x()
              * secondMoment(p1Tri.x(), p2Tri.x(), p3Tri.x())
              / S(12.0);
            integrals.y2 +=
                crossProdTri.y()
              * secondMoment(p1Tri.y(), p2Tri.y(), p3Tri.y())
              / S(12.0);
            integrals.z2 +=
                crossProdTri.z()
              * secondMoment(p1Tri.z(), p2Tri.z(), p3Tri.z())
              / S(12.0);

            Vector triangleCentroid = (p1Tri + p2Tri + p3Tri) / S(3.0);

            // Volume contribution for this sub-triangle
            integrals.volume += dot(triangleCentroid, crossProdTri) / S(2.0);

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

    return integrals;
}

// ************************* Distance Property Methods ************************

void Face::calculateDistanceProperties(std::span<const Vector> cellCentroids)
{
    dPf_ = centroid_ - cellCentroids[ownerCell_];
    dPfMag_ = dPf_.magnitude();

    // Calculate dNf only for internal faces
    if (!isBoundary())
    {
        size_t N = neighborCell_.value();

        Vector dNfVec = centroid_ - cellCentroids[N];
        dNf_ = dNfVec;
        dNfMag_ = dNfVec.magnitude();
    }
    distancePropertiesCalculated_ = true;
}

// ***************************** Output Methods *******************************

std::ostream& operator<<(std::ostream& os, const Face& f)
{
    os  << "Face(ID: " << f.idx() << ", Nodes: [";

    const auto& nodes = f.nodeIndices();
    for (size_t nodeIdx = 0; nodeIdx < nodes.size(); ++nodeIdx)
    {
        os  << nodes[nodeIdx]
            << (nodeIdx == nodes.size() - 1 ? "" : ", ");
    }

    os  <<  "], Owner: " << f.ownerCell() << ", Neighbor: "
        <<  (
                f.isBoundary() ? "Boundary"
              : std::to_string(f.neighborCell().value())
            );

    if (f.geometricPropertiesCalculated())
    {
        // save the current format
        auto flags = os.flags();
        auto prec = os.precision();

        // change format for geometric properties
        os  << std::fixed << std::setprecision(6);

        // Output geometric properties
        os  << ", Centroid: " << f.centroid()
            << ", Area: "   << f.projectedArea()
            << ", Normal: " << f.normal();

        // restore original format
        os.flags(flags);
        os.precision(prec);
    }
    else
    {
        os  << ", Geometry: N/A";
    }

    if (f.distancePropertiesCalculated())
    {
        // save the current format
        auto flags = os.flags();
        auto prec = os.precision();

        // change format for distance properties
        os  << std::fixed << std::setprecision(6);

        // Output distance properties
        os  << ", dPfMag: " << f.dPfMag();

        if (f.dNfMag().has_value())
        {
            os  << ", dNfMag: " << f.dNfMag().value();
        }

        // restore original format
        os.flags(flags);
        os.precision(prec);
    }

    os  << ')';

    return os;
}
