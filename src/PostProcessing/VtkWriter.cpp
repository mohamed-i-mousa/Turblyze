/******************************************************************************
 * @file VtkWriter.cpp
 * @brief Implementation of VTK output writer for post-processing
 *****************************************************************************/

// VTK includes
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#include <vtkPyramid.h>
#include <vtkIdList.h>

// STL includes
#include <stdexcept>
#include <set>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// Header includes
#include "VtkWriter.hpp"
#include "Cell.hpp"


namespace VtkWriter
{

// ****************** Private Helper Methods: Field Utilities *****************

ScalarField computeMagnitude
(
    const VectorField& field,
    const std::string& fieldName
)
{
    ScalarField magnitude(fieldName, field.size());
    for (size_t idx = 0; idx < field.size(); ++idx)
    {
        magnitude[idx] = field[idx].magnitude();
    }
    return magnitude;
}

StrainTensor computeStrainTensor
(
    const Vector& gradUx,
    const Vector& gradUy,
    const Vector& gradUz
)
{
    StrainTensor S;

    // Extract velocity gradient tensor components
    Scalar dUxdx = gradUx.x();
    Scalar dUxdy = gradUx.y();
    Scalar dUxdz = gradUx.z();

    Scalar dUydx = gradUy.x();
    Scalar dUydy = gradUy.y();
    Scalar dUydz = gradUy.z();

    Scalar dUzdx = gradUz.x();
    Scalar dUzdy = gradUz.y();
    Scalar dUzdz = gradUz.z();

    // Symmetric strain rate tensor: S_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
    S.S_11 = dUxdx;
    S.S_22 = dUydy;
    S.S_33 = dUzdz;
    S.S_12 = 0.5 * (dUxdy + dUydx);
    S.S_13 = 0.5 * (dUxdz + dUzdx);
    S.S_23 = 0.5 * (dUydz + dUzdy);

    return S;
}

// *************** Private Helper Methods: Geometric Utilities ****************

Vector computeFaceCentroid
(
    const std::vector<size_t>& faceNodes,
    const std::vector<Vector>& allNodes
)
{
    Vector centroid(0.0, 0.0, 0.0);
    for (size_t nodeIdx : faceNodes)
    {
        centroid += allNodes[nodeIdx];
    }
    return centroid / S(faceNodes.size());
}

Vector computeTriangleNormal
(
    const std::vector<size_t>& triNodes,
    const std::vector<Vector>& allNodes
)
{
    Vector v0 = allNodes[triNodes[0]];
    Vector v1 = allNodes[triNodes[1]];
    Vector v2 = allNodes[triNodes[2]];
    Vector edge1 = v1 - v0;
    Vector edge2 = v2 - v0;
    return cross(edge1, edge2);
}

Vector computeQuadNormal
(
    const std::vector<size_t>& quadNodes,
    const std::vector<Vector>& allNodes
)
{
    Vector v0 = allNodes[quadNodes[0]];
    Vector v1 = allNodes[quadNodes[1]];
    Vector v2 = allNodes[quadNodes[2]];
    Vector edge1 = v1 - v0;
    Vector edge2 = v2 - v0;
    return cross(edge1, edge2);
}

// ************ Private Helper Methods: VTK Cell Node Ordering ****************

std::vector<vtkIdType> orderHexahedronNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::set<size_t>& uniqueNodes,
    const std::vector<Vector>& allNodes
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
            return orderedNodes; // Not a standard hex
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
    size_t bottomFaceIdx = 0, topFaceIdx = 0;

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
    Vector bottomNormal = computeQuadNormal(bottomFace, allNodes);
    Vector topNormal = computeQuadNormal(topFace, allNodes);
    Vector bottomCentroid = computeFaceCentroid(bottomFace, allNodes);
    Vector topCentroid = computeFaceCentroid(topFace, allNodes);

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
    Vector testNormal = computeQuadNormal(orderedBottomNodes, allNodes);

    Vector outwardDir =
        computeFaceCentroid(orderedBottomNodes, allNodes) - centroid;

    if (dot(testNormal, outwardDir) < 0)
    {
        std::reverse(orderedBottomNodes.begin(), orderedBottomNodes.end());
    }

    // Build edge connectivity from side faces to match top to bottom
    // Each side face has 4 nodes forming a quad: 2 from bottom, 2 from top
    // We need to find which bottom node connects to which top node via
    // vertical edges
    std::map<size_t, size_t> nodeConnections;

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

        if (nodeConnections.find(bottomNode) != nodeConnections.end())
        {
            // Direct one-to-one mapping from bottom to top node
            size_t topNode = nodeConnections[bottomNode];

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
        }
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
    const std::set<size_t>& uniqueNodes,
    const std::vector<Vector>& allNodes
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
    Vector tri0Centroid = computeFaceCentroid(triangularFaces[0], allNodes);
    Vector tri1Centroid = computeFaceCentroid(triangularFaces[1], allNodes);

    Vector tri0Normal = computeTriangleNormal(triangularFaces[0], allNodes);
    Vector tri1Normal = computeTriangleNormal(triangularFaces[1], allNodes);

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
    std::map<size_t, size_t> bottomToTopMap;

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
    Vector testNormal = computeTriangleNormal(orderedBottomNodes, allNodes);

    Vector outwardDir =
        computeFaceCentroid(orderedBottomNodes, allNodes) - centroid;
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
    const std::set<size_t>& uniqueNodes,
    const std::vector<Vector>& allNodes
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
        computeQuadNormal(orderedBase, allNodes);

    Vector outwardDir =
        computeFaceCentroid(orderedBase, allNodes) - centroid;

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

// ******************* Public API: Derived Field Computation *******************

ScalarField computeVelocityMagnitude(const VectorField& velocity)
{
    return computeMagnitude(velocity, "velocityMagnitude");
}

ScalarField computeVorticityMagnitude(const VectorField& vorticity)
{
    return computeMagnitude(vorticity, "vorticityMagnitude");
}

ScalarField computeQCriterion
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
)
{
    ScalarField qCriterion("QCriterion", gradUx.size());

    for (size_t i = 0; i < gradUx.size(); ++i)
    {
        // Compute strain rate tensor using helper function
        StrainTensor S = computeStrainTensor(gradUx[i], gradUy[i], gradUz[i]);

        // Extract velocity gradient tensor components for rotation tensor
        Scalar dUxdy = gradUx[i].y();
        Scalar dUxdz = gradUx[i].z();
        Scalar dUydx = gradUy[i].x();
        Scalar dUydz = gradUy[i].z();
        Scalar dUzdx = gradUz[i].x();
        Scalar dUzdy = gradUz[i].y();

        // Anti-symmetric rotation rate tensor
        // (Omega_ij = 0.5 * (du_i/dx_j - du_j/dx_i))
        Scalar Omega_12 = 0.5 * (dUxdy - dUydx);
        Scalar Omega_13 = 0.5 * (dUxdz - dUzdx);
        Scalar Omega_23 = 0.5 * (dUydz - dUzdy);

        // Q-criterion = 0.5 * (||Omega||^2 - ||S||^2)
        Scalar omegaNormSq =
            2.0 * (Omega_12*Omega_12 + Omega_13*Omega_13 + Omega_23*Omega_23);

        qCriterion[i] = 0.5 * (omegaNormSq - S.normSquared());
    }

    return qCriterion;
}

ScalarField computeStrainRateMagnitude
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
)
{
    ScalarField strainRateMag("strainRateMagnitude", gradUx.size());

    for (size_t idx = 0; idx < gradUx.size(); ++idx)
    {
        // Compute strain rate tensor using helper function
        StrainTensor S = computeStrainTensor(gradUx[idx], gradUy[idx], gradUz[idx]);

        // Strain rate magnitude = sqrt(2 * S_ij * S_ij)
        strainRateMag[idx] = std::sqrt(2.0 * S.normSquared());
    }

    return strainRateMag;
}

// *********************** Public API: VTK File Export ************************

void writeVtkUnstructuredGrid
(
    const std::string& filename,
    const std::vector<Vector>& allNodes,
    const std::vector<Cell>& allCells,
    const std::vector<Face>& allFaces,
    const std::map<std::string,
    const ScalarField*>& scalarCellFields,
    const std::map<std::string,
    const VectorField*>& vectorCellFields
)
{
    // Create vtkUnstructuredGrid
    vtkSmartPointer<vtkUnstructuredGrid>
    unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Create vtkPoints from nodes
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(allNodes.size());

    for (size_t i = 0; i < allNodes.size(); ++i)
    {
        points->SetPoint(i, allNodes[i].x(), allNodes[i].y(), allNodes[i].z());
    }

    unstructuredGrid->SetPoints(points);

    // Add cells to unstructured grid with proper topology preservation
    // Build a node connectivity graph to determine proper VTK node ordering

    // Cell type statistics for diagnostics
    size_t numTets = 0, numHexes = 0, numWedges = 0,
           numPyramids = 0, numConvex = 0;
    size_t numWedgeOrderingFailures = 0;

    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        const Cell& cell = allCells[cellIdx];
        const auto& faceIndices = cell.faceIndices();

        // Collect all unique nodes and build connectivity
        std::set<size_t> uniqueNodeSet;
        std::vector<std::vector<size_t>> faceNodeLists;

        for (size_t faceIdx : faceIndices)
        {
            const Face& face = allFaces[faceIdx];
            const auto& nodeIndices = face.nodeIndices();
            uniqueNodeSet.insert(nodeIndices.begin(), nodeIndices.end());
            faceNodeLists.push_back(nodeIndices);
        }

        size_t numNodes = uniqueNodeSet.size();
        size_t numFaces = faceIndices.size();

        // Determine cell type
        int cellType = VTK_CONVEX_POINT_SET; // Fallback
        std::vector<vtkIdType> orderedNodes;

        if (numNodes == 4 && numFaces == 4)
        {
            // Tetrahedron - VTK ordering:
            // any vertex, then 3 neighbors forming base
            cellType = VTK_TETRA;
            numTets++;
            orderedNodes.reserve(4);
            for (size_t nodeIdx : uniqueNodeSet)
            {
                orderedNodes.push_back(static_cast<vtkIdType>(nodeIdx));
            }
        }
        else if (numNodes == 8 && numFaces == 6)
        {
            // Hexahedron - need to find proper ordering
            cellType = VTK_HEXAHEDRON;
            numHexes++;
            orderedNodes =
                orderHexahedronNodes(faceNodeLists, uniqueNodeSet, allNodes);
        }
        else if (numNodes == 6 && numFaces == 5)
        {
            // Wedge (prism)
            cellType = VTK_WEDGE;
            numWedges++;
            orderedNodes =
                orderWedgeNodes(faceNodeLists, uniqueNodeSet, allNodes);

            if (orderedNodes.empty())
            {
                numWedgeOrderingFailures++;
                // Fallback to convex point set
                cellType = VTK_CONVEX_POINT_SET;
                for (size_t nodeIdx : uniqueNodeSet)
                {
                    orderedNodes.push_back(static_cast<vtkIdType>(nodeIdx));
                }
            }
        }
        else if (numNodes == 5 && numFaces == 5)
        {
            // Pyramid
            cellType = VTK_PYRAMID;
            numPyramids++;
            orderedNodes =
                orderPyramidNodes(faceNodeLists, uniqueNodeSet, allNodes);
        }
        else
        {
            // Use convex point set for unknown topology
            numConvex++;
            orderedNodes.reserve(numNodes);
            for (size_t nodeId : uniqueNodeSet)
            {
                orderedNodes.push_back(static_cast<vtkIdType>(nodeId));
            }
        }

        // Insert cell with ordered nodes
        if (orderedNodes.empty())
        {
            // Fallback: use arbitrary ordering
            vtkSmartPointer<vtkIdList> cellPointIds =
                vtkSmartPointer<vtkIdList>::New();

            for (size_t nodeIdx : uniqueNodeSet)
            {
                cellPointIds->InsertNextId(static_cast<vtkIdType>(nodeIdx));
            }
            unstructuredGrid->InsertNextCell
            (
                VTK_CONVEX_POINT_SET,
                cellPointIds
            );
        }
        else
        {
            unstructuredGrid->InsertNextCell
            (
                cellType,
                static_cast<vtkIdType>(orderedNodes.size()),
                orderedNodes.data()
            );
        }
    }

    // Add scalar cell fields
    for (const auto& [fieldName, scalarField] : scalarCellFields)
    {
        if (!scalarField) continue;

        vtkSmartPointer<vtkDoubleArray> 
        dataArray = vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(allCells.size());

        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            if (cellIdx < scalarField->size())
            {
                dataArray->SetValue(cellIdx, (*scalarField)[cellIdx]);
            }
            else
            {
                dataArray->SetValue(cellIdx, 0.0);
            }
        }

        unstructuredGrid->GetCellData()->AddArray(dataArray);
    }

    // Add vector cell fields
    for (const auto& [fieldName, vectorField] : vectorCellFields)
    {
        if (!vectorField) continue;

        vtkSmartPointer<vtkDoubleArray> dataArray =
            vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(3);
        dataArray->SetNumberOfTuples(allCells.size());

        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            if (cellIdx < vectorField->size())
            {
                const Vector& vec = (*vectorField)[cellIdx];
                double vecData[3] = {vec.x(), vec.y(), vec.z()};
                dataArray->SetTuple(cellIdx, vecData);
            }
            else
            {
                double zeroVec[3] = {0.0, 0.0, 0.0};
                dataArray->SetTuple(cellIdx, zeroVec);
            }
        }

        unstructuredGrid->GetCellData()->AddArray(dataArray);
    }

    // Write the unstructured grid to file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    writer->SetFileName(filename.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->SetDataModeToBinary();
    writer->SetCompressorTypeToZLib();
    writer->Write();

    std::cout
        << "VTK UnstructuredGrid file written successfully: "
        << filename << std::endl;

    std::cout
        << "  - Number of points: " << allNodes.size() << std::endl;

    std::cout
        << "  - Number of cells: " << allCells.size() << std::endl;

    std::cout
        << "  - Cell types:" << std::endl;

    if (numTets > 0)
    {
        std::cout
            << "    - Tetrahedra: " << numTets << std::endl;
    }

    if (numHexes > 0)
    {
        std::cout
            << "    - Hexahedra: " << numHexes << std::endl;
    }

    if (numWedges > 0)
    {
        std::cout
            << "    - Wedges (prisms): " << numWedges << std::endl;
    }

    if (numPyramids > 0)
    {
        std::cout
            << "    - Pyramids: " << numPyramids << std::endl;
    }

    if (numConvex > 0)
    {
        std::cout
            << "    - Convex point sets: " << numConvex << std::endl;
    }

    if (numWedgeOrderingFailures > 0)
    {
        std::cout
            << "  - WARNING: " << numWedgeOrderingFailures
            << " wedge cells failed proper ordering"
            << std::endl;
    }

    std::cout
        << "  - Number of scalar fields: " << scalarCellFields.size()
        << std::endl;

    std::cout
        << "  - Number of vector fields: " << vectorCellFields.size()
        << std::endl;
}

// *********************** Public API: PVD Time Series ************************

void writePVDTimeSeriesHeader
(
    const std::string& pvdFilename
)
{
    std::ofstream pvdFile(pvdFilename);

    if (!pvdFile.is_open())
    {
        throw   std::runtime_error
                (
                    "Failed to open PVD file: " + pvdFilename
                );
    }

    pvdFile << "<?xml version=\"1.0\"?>\n";
    pvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" "
            << "byte_order=\"LittleEndian\">\n";
    pvdFile << "  <Collection>\n";
    pvdFile << "  </Collection>\n";
    pvdFile << "</VTKFile>\n";

    pvdFile.close();

    std::cout
        << "PVD time series header written: "
        << pvdFilename << std::endl;
}

void appendPVDTimeStep
(
    const std::string& pvdFilename,
    const std::string& vtuFilename,
    Scalar timeValue
)
{
    // Read existing PVD file
    std::ifstream pvdFileIn(pvdFilename);
    if (!pvdFileIn.is_open())
    {
        throw   std::runtime_error
                (
                    "Failed to open PVD file for reading: " + pvdFilename
                );
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(pvdFileIn, line))
    {
        lines.push_back(line);
    }
    pvdFileIn.close();

    // Find the </Collection> line and insert before it
    std::ofstream pvdFileOut(pvdFilename);
    if (!pvdFileOut.is_open())
    {
        throw   std::runtime_error
                (
                    "Failed to open PVD file for writing: " + pvdFilename
                );
    }

    for (const auto& existingLine : lines)
    {
        if (existingLine.find("</Collection>") != std::string::npos)
        {
            // Insert the new timestep before closing collection
            pvdFileOut << "    <DataSet timestep=\"" << timeValue
                      << "\" file=\"" << vtuFilename << "\"/>\n";
        }
        pvdFileOut << existingLine << "\n";
    }

    pvdFileOut.close();

    std::cout
        << "Added timestep " << timeValue << " to PVD file"
        << std::endl;
}

// *********************** Public API: Geometry Export ************************

void writeCellGeometryData
(
    const std::string& filename,
    const std::vector<Cell>& allCells
)
{
    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        throw   std::runtime_error
                (
                    "Failed to open file: " + filename
                );
    }

    // Header line
    outFile << "# cellIndex volume centroid_x centroid_y centroid_z\n";
    outFile << std::scientific << std::setprecision(12);

    for (size_t i = 0; i < allCells.size(); ++i)
    {
        const Cell& cell = allCells[i];
        const Vector& c = cell.centroid();

        outFile << i << " "
                << cell.volume() << " "
                << c.x() << " "
                << c.y() << " "
                << c.z() << "\n";
    }

    outFile.close();
    std::cout
        << "Cell geometry data written to: " << filename << std::endl;

    std::cout
        << "  - Number of cells: " << allCells.size() << std::endl;
}
} // namespace VtkWriter