/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Face.cpp
 * @brief Implementation of face geometric properties and operations
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Face.h"

// Standard library headers
#include <cmath>

// Project headers
#include "Cell.h"
#include "ErrorHandler.h"

// ***************************** Internal Helpers *****************************

namespace
{

// Symmetric second-moment polynomial for triangle integration
// Evaluates a² + b² + c² + ab + ac + bc
// ∫∫_triangle x² dA = (area / 6) × secondMoment(x₁, x₂, x₃)
[[nodiscard]] Scalar secondMoment
(
    Scalar a,
    Scalar b,
    Scalar c
) noexcept
{
    return a*a + b*b + c*c + a*b + a*c + b*c;
}

} // namespace

// ************************ Geometric Property Methods ************************

FaceIntegrals Face::geometricProperties
(
    const NodeList& allNodes
)
{
    const Count numNodes = nodeIndices_.size();

    for (Index nodeIdx : nodeIndices_)
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

        const Vector vecA = p2 - p1;
        const Vector vecB = p3 - p1;
        const Vector crossProd = cross(vecB, vecA);
        const Scalar crossProdMag = magnitude(crossProd);

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
    }
    // CASE 2: Face is "Polygon" (numNodes > 3)
    else
    {
        Vector faceCenter;

        for (Index nodeIdx : nodeIndices_)
        {
            faceCenter += allNodes[nodeIdx];
        }

        faceCenter /= S(numNodes);

        Vector weightedCentroidSum{};
        Vector normalSum{};
        Scalar weightedAreaSum = S(0.0);

        for (Index nodeIdx = 0; nodeIdx < numNodes; ++nodeIdx)
        {
            const Vector& pCurr = allNodes[nodeIndices_[nodeIdx]];
            const Vector& pNext =
                allNodes[nodeIndices_[(nodeIdx + 1) % numNodes]];

            const Vector& p1Tri = faceCenter;
            const Vector& p2Tri = pCurr;
            const Vector& p3Tri = pNext;

            const Vector vecATri = p2Tri - p1Tri;
            const Vector vecBTri = p3Tri - p1Tri;
            const Vector crossProdTri = cross(vecBTri, vecATri);
            const Scalar triangleArea = S(0.5) * magnitude(crossProdTri);
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

            const Vector triangleCentroid =
                (p1Tri + p2Tri + p3Tri) / S(3.0);

            // Volume contribution for this sub-triangle
            integrals.volume += dot(triangleCentroid, crossProdTri) / S(2.0);

            weightedAreaSum += triangleArea;
            weightedCentroidSum += triangleCentroid * triangleArea;
        }

        projectedArea_ = magnitude(normalSum) / S(2.0);

        contactArea_ = weightedAreaSum;

        // Centroid uses contact area weighting (sum of triangle areas)
        centroid_ = weightedCentroidSum / (weightedAreaSum + vSmallValue);

        normal_ = normalized(normalSum);
    }

    return integrals;
}
