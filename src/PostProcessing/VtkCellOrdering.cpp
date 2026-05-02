/******************************************************************************
 * @file VtkCellOrdering.cpp
 * @brief Implementation of VTK cell node-ordering helpers
 *****************************************************************************/

#include "VtkCellOrdering.h"

#include <algorithm>
#include <set>
#include <unordered_map>

#include "Scalar.h"


namespace VTK
{

// *************** Internal Helper Methods: Geometric Utilities ****************

Vector faceCentroid
(
    std::span<const size_t> faceNodes,
    std::span<const Vector> allNodes
)
{
    Vector centroid(0.0, 0.0, 0.0);
    for (size_t nodeIdx : faceNodes)
    {
        centroid += allNodes[nodeIdx];
    }
    return centroid / S(faceNodes.size());
}

Vector faceNormal
(
    std::span<const size_t> faceNodes,
    std::span<const Vector> allNodes
)
{
    Vector v0 = allNodes[faceNodes[0]];
    Vector v1 = allNodes[faceNodes[1]];
    Vector v2 = allNodes[faceNodes[2]];
    Vector edge1 = v1 - v0;
    Vector edge2 = v2 - v0;
    return cross(edge1, edge2);
}

// ************ Internal Helper Methods: VTK Cell Node Ordering ****************

std::vector<vtkIdType> orderHexahedronNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::unordered_set<size_t>& uniqueNodes,
    std::span<const Vector> allNodes
)
{
    std::vector<vtkIdType> orderedNodes;

    if (uniqueNodes.size() != 8 || faceNodeLists.size() != 6)
    {
        return orderedNodes; // Invalid hex
    }

    // All faces should be quads
    for (const auto& faceNodes : faceNodeLists)
    {
        if (faceNodes.size() != 4)
        {
            return orderedNodes; // Not a hex
        }
    }

    // Compute cell centroid
    Vector centroid(0.0, 0.0, 0.0);
    for (size_t nodeId : uniqueNodes)
    {
        centroid += allNodes[nodeId];
    }
    centroid /= 8.0;

    // Find two opposite faces (bottom and top)
    // Opposite faces share no common nodes
    std::vector<size_t> bottomFace, topFace;
    size_t bottomFaceIdx = 0;
    size_t topFaceIdx = 0;

    for (size_t i = 0; i < faceNodeLists.size(); ++i)
    {
        for (size_t j = i + 1; j < faceNodeLists.size(); ++j)
        {
            const auto& face1 = faceNodeLists[i];
            const auto& face2 = faceNodeLists[j];

            // Check if faces share no nodes (opposite faces)
            bool shareNode = false;
            for (size_t n1 : face1)
            {
                for (size_t n2 : face2)
                {
                    if (n1 == n2)
                    {
                        shareNode = true;
                        break;
                    }
                }
                if (shareNode) break;
            }

            if (!shareNode)
            {
                bottomFace = face1;
                topFace = face2;
                bottomFaceIdx = i;
                topFaceIdx = j;
                break;
            }
        }
        if (!bottomFace.empty()) break;
    }

    if (bottomFace.empty() || topFace.empty())
    {
        return orderedNodes; // Failed to find opposite faces
    }

    // Determine face orientation using face normals
    Vector bottomNormal = faceNormal(bottomFace, allNodes);
    Vector topNormal = faceNormal(topFace, allNodes);
    Vector bottomCentroid = faceCentroid(bottomFace, allNodes);
    Vector topCentroid = faceCentroid(topFace, allNodes);

    Vector bottomToCenter = centroid - bottomCentroid;
    Vector topToCenter = centroid - topCentroid;

    // Orient so bottom face is the one with normal pointing inward
    if (dot(bottomNormal, bottomToCenter) < dot(topNormal, topToCenter))
    {
        std::swap(bottomFace, topFace);
    }

    // Order bottom face nodes consistently
    std::vector<size_t> orderedBottomNodes = bottomFace;

    // Check winding and reverse if needed
    Vector testNormal = faceNormal(orderedBottomNodes, allNodes);

    Vector outwardDir =
        faceCentroid(orderedBottomNodes, allNodes) - centroid;

    if (dot(testNormal, outwardDir) < 0)
    {
        std::reverse(orderedBottomNodes.begin(), orderedBottomNodes.end());
    }

    // Build edge connectivity from side faces to match top to bottom
    // Each side face has 4 nodes forming a quad: 2 from bottom, 2 from top
    // We need to find which bottom node connects to which top node via
    // vertical edges
    std::unordered_map<size_t, size_t> nodeConnections;

    for (size_t i = 0; i < faceNodeLists.size(); ++i)
    {
        if (i == bottomFaceIdx || i == topFaceIdx) continue; // Skip top/bottom

        const auto& sideFace = faceNodeLists[i];
        if (sideFace.size() != 4) continue;

        // For each edge of the side face, check if it connects
        // a bottom node to a top node
        for (size_t e = 0; e < 4; ++e)
        {
            size_t n0 = sideFace[e];
            size_t n1 = sideFace[(e + 1) % 4];

            bool n0_bottom =
                std::find
                (bottomFace.begin(), bottomFace.end(), n0)
              != bottomFace.end();

            bool n1_bottom =
                std::find
                (bottomFace.begin(), bottomFace.end(), n1)
              != bottomFace.end();

            bool n0_top =
                std::find(topFace.begin(), topFace.end(), n0) != topFace.end();

            bool n1_top =
                std::find(topFace.begin(), topFace.end(), n1) != topFace.end();

            // If the edge connects a bottom node to a top node (vertical edge)
            if (n0_bottom && n1_top)
            {
                nodeConnections[n0] = n1;
            }
            else if (n1_bottom && n0_top)
            {
                nodeConnections[n1] = n0;
            }
        }
    }

    // Match top nodes to bottom nodes
    std::vector<size_t> orderedTopNodes(4);
    for (size_t i = 0; i < 4; ++i)
    {
        size_t bottomNode = orderedBottomNodes[i];

        auto it = nodeConnections.find(bottomNode);
        if (it == nodeConnections.end())
        {
            // Failed to find mapping — return empty to trigger fallback
            return std::vector<vtkIdType>();
        }

        size_t topNode = it->second;

        // Verify it's actually in the top face
        if
        (
            std::find(topFace.begin(),
            topFace.end(),
            topNode) != topFace.end()
        )
        {
            orderedTopNodes[i] = topNode;
        }
        else
        {
            return std::vector<vtkIdType>();
        }
    }

    // Verify all top nodes were mapped correctly
    std::set<size_t> mappedTopNodes
    (
        orderedTopNodes.begin(),
        orderedTopNodes.end()
    );

    std::set<size_t> expectedTopNodes(topFace.begin(), topFace.end());

    if (mappedTopNodes != expectedTopNodes)
    {
        return std::vector<vtkIdType>();
    }

    // Build final VTK hex ordering
    orderedNodes.reserve(8);
    for (size_t nodeIdx : orderedBottomNodes)
    {
        orderedNodes.push_back(static_cast<vtkIdType>(nodeIdx));
    }
    for (size_t nodeIdx : orderedTopNodes)
    {
        orderedNodes.push_back(static_cast<vtkIdType>(nodeIdx));
    }

    return orderedNodes;
}

std::vector<vtkIdType> orderWedgeNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::unordered_set<size_t>& uniqueNodes,
    std::span<const Vector> allNodes
)
{
    std::vector<vtkIdType> orderedNodes;

    if (uniqueNodes.size() != 6 || faceNodeLists.size() != 5)
    {
        return orderedNodes;
    }

    // Identify triangular faces (3 nodes) and quadrilateral faces (4 nodes)
    std::vector<std::vector<size_t>> triangularFaces;
    std::vector<std::vector<size_t>> quadFaces;

    for (const auto& faceNodes : faceNodeLists)
    {
        if (faceNodes.size() == 3)
        {
            triangularFaces.push_back(faceNodes);
        }
        else if (faceNodes.size() == 4)
        {
            quadFaces.push_back(faceNodes);
        }
    }

    // Prism must have exactly 2 triangular faces and 3 quad faces
    if (triangularFaces.size() != 2 || quadFaces.size() != 3)
    {
        return orderedNodes;
    }

    // Determine which triangular face is "bottom" and which is "top"
    // Use centroid as reference - doesn't assume Z-alignment
    Vector centroid(0.0, 0.0, 0.0);
    for (size_t nodeId : uniqueNodes)
    {
        centroid += allNodes[nodeId];
    }
    centroid /= 6.0;

    // Compute face centroids and normals for triangular faces
    Vector tri0Centroid = faceCentroid(triangularFaces[0], allNodes);
    Vector tri1Centroid = faceCentroid(triangularFaces[1], allNodes);

    Vector tri0Normal = faceNormal(triangularFaces[0], allNodes);
    Vector tri1Normal = faceNormal(triangularFaces[1], allNodes);

    // Check which triangle is more "inward" pointing relative to centroid
    Vector tri0ToCenter = centroid - tri0Centroid;
    Vector tri1ToCenter = centroid - tri1Centroid;

    Scalar dot0 = dot(tri0Normal, tri0ToCenter);
    Scalar dot1 = dot(tri1Normal, tri1ToCenter);

    // The face whose normal points toward the centroid is the "bottom"
    std::vector<size_t> bottomTri, topTri;
    if (dot0 > dot1)
    {
        bottomTri = triangularFaces[0];
        topTri = triangularFaces[1];
    }
    else
    {
        bottomTri = triangularFaces[1];
        topTri = triangularFaces[0];
    }

    // Build edge connectivity from quad faces to establish correspondence
    // In a prism, each quad face has 4 nodes forming a cycle
    // Two nodes are from bottom triangle (forming an edge), two from top
    // triangle (forming an edge)
    // The edges are connected by the quad face
    std::unordered_map<size_t, size_t> bottomToTopMap;

    for (const auto& quadNodes : quadFaces)
    {
        if (quadNodes.size() != 4) continue;

        // Classify nodes as bottom or top
        std::vector<size_t> quadBottomNodes, quadTopNodes;

        for (size_t nid : quadNodes)
        {
            bool inBottom =
                std::find
                (bottomTri.begin(), bottomTri.end(), nid)
              != bottomTri.end();
            if (inBottom)
            {
                quadBottomNodes.push_back(nid);
            }
            else
            {
                quadTopNodes.push_back(nid);
            }
        }

        if (quadBottomNodes.size() != 2 || quadTopNodes.size() != 2) continue;

        // Find which nodes in the quad are adjacent (form edges)
        // In the quad [n0, n1, n2, n3], edges are:
        // n0-n1, n1-n2, n2-n3, n3-n0

        // Check all 4 edges of the quad to find which bottom and top
        // nodes are connected
        for (size_t i = 0; i < 4; ++i)
        {
            size_t n0 = quadNodes[i];
            size_t n1 = quadNodes[(i + 1) % 4];

            bool n0_bottom =
                std::find
                (bottomTri.begin(), bottomTri.end(), n0)
             != bottomTri.end();

            bool n1_bottom =
                std::find(bottomTri.begin(), bottomTri.end(), n1)
             != bottomTri.end();

            // If this edge connects a bottom node to a top node,
            // record the connection
            if (n0_bottom && !n1_bottom)
            {
                bottomToTopMap[n0] = n1;
            }
            else if (!n0_bottom && n1_bottom)
            {
                bottomToTopMap[n1] = n0;
            }
        }
    }

    // Order bottom triangle nodes consistently
    // (counter-clockwise when viewed from outside)
    std::vector<size_t> orderedBottomNodes = bottomTri;

    // Ensure proper winding by checking normal direction
    Vector testNormal = faceNormal(orderedBottomNodes, allNodes);

    Vector outwardDir =
        faceCentroid(orderedBottomNodes, allNodes) - centroid;
    if (dot(testNormal, outwardDir) < 0)
    {
        // Reverse winding to make it outward-pointing
        std::reverse(orderedBottomNodes.begin(), orderedBottomNodes.end());
    }

    // Match top nodes to bottom nodes using the mapping we built
    std::vector<size_t> orderedTopNodes(3);

    for (size_t i = 0; i < 3; ++i)
    {
        size_t bottomNode = orderedBottomNodes[i];

        // Look up the corresponding top node
        if (bottomToTopMap.find(bottomNode) != bottomToTopMap.end())
        {
            orderedTopNodes[i] = bottomToTopMap[bottomNode];
        }
        else
        {
            // Failed to find mapping - return empty to trigger fallback
            return std::vector<vtkIdType>();
        }
    }

    // Verify all top nodes were mapped and are valid
    std::set<size_t> mappedTopNodes
    (
        orderedTopNodes.begin(),
        orderedTopNodes.end()
    );

    std::set<size_t> expectedTopNodes(topTri.begin(), topTri.end());

    if (mappedTopNodes != expectedTopNodes)
    {
        // Mapping is incomplete or incorrect. Return empty to trigger fallback
        return std::vector<vtkIdType>();
    }

    // Build final VTK wedge ordering
    orderedNodes.reserve(6);
    for (size_t nodeId : orderedBottomNodes)
    {
        orderedNodes.push_back(static_cast<vtkIdType>(nodeId));
    }
    for (size_t nodeId : orderedTopNodes)
    {
        orderedNodes.push_back(static_cast<vtkIdType>(nodeId));
    }

    return orderedNodes;
}

std::vector<vtkIdType> orderPyramidNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::unordered_set<size_t>& uniqueNodes,
    std::span<const Vector> allNodes
)
{
    std::vector<vtkIdType> orderedNodes;

    if (uniqueNodes.size() != 5 || faceNodeLists.size() != 5)
    {
        return orderedNodes;
    }

    // Identify the quad face (base) and triangular faces
    // A pyramid has 1 quad face and 4 triangular faces
    std::vector<size_t> baseFace;

    for (const auto& faceNodes : faceNodeLists)
    {
        if (faceNodes.size() == 4)
        {
            baseFace = faceNodes;
            break;
        }
    }

    if (baseFace.size() != 4)
    {
        return orderedNodes;
    }

    // The apex is the node NOT in the quad base face
    size_t apexNode = 0;
    for (size_t nodeId : uniqueNodes)
    {
        if
        (
            std::find
            (
                baseFace.begin(),
                baseFace.end(),
                nodeId
            )
         == baseFace.end()
        )
        {
            apexNode = nodeId;
            break;
        }
    }

    // Compute cell centroid
    Vector centroid(0.0, 0.0, 0.0);
    for (size_t nodeId : uniqueNodes)
    {
        centroid += allNodes[nodeId];
    }
    centroid /= 5.0;

    // Ensure proper winding of base quad
    // (outward-pointing normal when viewed from outside)
    std::vector<size_t> orderedBase = baseFace;
    Vector testNormal =
        faceNormal(orderedBase, allNodes);

    Vector outwardDir =
        faceCentroid(orderedBase, allNodes) - centroid;

    if (dot(testNormal, outwardDir) < 0)
    {
        std::reverse(orderedBase.begin(), orderedBase.end());
    }

    // Build ordered list: 4 base nodes + apex
    orderedNodes.reserve(5);
    for (size_t nodeId : orderedBase)
    {
        orderedNodes.push_back(
            static_cast<vtkIdType>(nodeId));
    }
    orderedNodes.push_back(
        static_cast<vtkIdType>(apexNode));

    return orderedNodes;
}

} // namespace VTK
