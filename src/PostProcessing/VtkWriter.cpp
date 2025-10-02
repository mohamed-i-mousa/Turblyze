/******************************************************************************
 * @file VtkWriter.cpp
 * @brief Implementation of VTK output writer for CFD results
 *****************************************************************************/

#include "VtkWriter.hpp"
#include "Cell.hpp"

// VTK includes
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#include <vtkPyramid.h>
#include <vtkIdList.h>

#include <stdexcept>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>

namespace VtkWriter
{

// Helper function to compute velocity magnitude
ScalarField computeVelocityMagnitude(const VectorField& velocity)
{
    ScalarField velocityMag("velocityMagnitude", velocity.size());

    for (size_t i = 0; i < velocity.size(); ++i)
    {
        velocityMag[i] = velocity[i].magnitude();
    }

    return velocityMag;
}

// Helper function to compute vorticity magnitude
ScalarField computeVorticityMagnitude(const VectorField& vorticity)
{
    ScalarField vorticityMag("vorticityMagnitude", vorticity.size());

    for (size_t i = 0; i < vorticity.size(); ++i)
    {
        vorticityMag[i] = vorticity[i].magnitude();
    }

    return vorticityMag;
}

// Helper function to compute Q-criterion
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
        // Velocity gradient tensor components
        Scalar dUxdx = gradUx[i].x();
        Scalar dUxdy = gradUx[i].y();
        Scalar dUxdz = gradUx[i].z();

        Scalar dUydx = gradUy[i].x();
        Scalar dUydy = gradUy[i].y();
        Scalar dUydz = gradUy[i].z();

        Scalar dUzdx = gradUz[i].x();
        Scalar dUzdy = gradUz[i].y();
        Scalar dUzdz = gradUz[i].z();

        // Symmetric strain rate tensor (S_ij = 0.5 * (du_i/dx_j + du_j/dx_i))
        Scalar S11 = dUxdx;
        Scalar S22 = dUydy;
        Scalar S33 = dUzdz;
        Scalar S12 = 0.5 * (dUxdy + dUydx);
        Scalar S13 = 0.5 * (dUxdz + dUzdx);
        Scalar S23 = 0.5 * (dUydz + dUzdy);

        // Anti-symmetric rotation rate tensor (Omega_ij = 0.5 * (du_i/dx_j - du_j/dx_i))
        Scalar Omega12 = 0.5 * (dUxdy - dUydx);
        Scalar Omega13 = 0.5 * (dUxdz - dUzdx);
        Scalar Omega23 = 0.5 * (dUydz - dUzdy);

        // Q-criterion = 0.5 * (||Omega||^2 - ||S||^2)
        Scalar omegaNormSq = 2.0 * (Omega12*Omega12 + Omega13*Omega13 + Omega23*Omega23);
        Scalar sNormSq = S11*S11 + S22*S22 + S33*S33 + 2.0*(S12*S12 + S13*S13 + S23*S23);

        qCriterion[i] = 0.5 * (omegaNormSq - sNormSq);
    }

    return qCriterion;
}

// Helper function to compute strain rate magnitude
ScalarField computeStrainRateMagnitude
(
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz
)
{
    ScalarField strainRateMag("strainRateMagnitude", gradUx.size());

    for (size_t i = 0; i < gradUx.size(); ++i)
    {
        // Velocity gradient tensor components
        Scalar dUxdx = gradUx[i].x();
        Scalar dUxdy = gradUx[i].y();
        Scalar dUxdz = gradUx[i].z();

        Scalar dUydx = gradUy[i].x();
        Scalar dUydy = gradUy[i].y();
        Scalar dUydz = gradUy[i].z();

        Scalar dUzdx = gradUz[i].x();
        Scalar dUzdy = gradUz[i].y();
        Scalar dUzdz = gradUz[i].z();

        // Symmetric strain rate tensor
        Scalar S11 = dUxdx;
        Scalar S22 = dUydy;
        Scalar S33 = dUzdz;
        Scalar S12 = 0.5 * (dUxdy + dUydx);
        Scalar S13 = 0.5 * (dUxdz + dUzdx);
        Scalar S23 = 0.5 * (dUydz + dUzdy);

        // Strain rate magnitude = sqrt(2 * S_ij * S_ij)
        Scalar sNormSq = S11*S11 + S22*S22 + S33*S33 + 2.0*(S12*S12 + S13*S13 + S23*S23);
        strainRateMag[i] = std::sqrt(2.0 * sNormSq);
    }

    return strainRateMag;
}

void writeVtkFile
(
    const std::string& filename,
    const std::vector<Vector>& allNodes,
    const std::vector<Face>& allFaces,
    const std::map<std::string,
    const ScalarField*>& scalarCellFields,
    const std::map<std::string,
    const VectorField*>& vectorCellFields
) 
{
    // Create a vtkPolyData object
    vtkSmartPointer<vtkPolyData> polyData = 
        vtkSmartPointer<vtkPolyData>::New();
    
    // Create vtkPoints from the allNodes vector
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    points->SetNumberOfPoints(allNodes.size());
    
    // Add all nodes to the points
    for (size_t i = 0; i < allNodes.size(); ++i) 
    {
        points->SetPoint(i, allNodes[i].x(), allNodes[i].y(), allNodes[i].z());
    }
    
    // Set the points to the polyData
    polyData->SetPoints(points);
    
    // Create vtkCellArray to store the polygons (faces)
    vtkSmartPointer<vtkCellArray> polygons =
        vtkSmartPointer<vtkCellArray>::New();
    
    // Keep track of owner cell for each face for mapping cell data to face data
    std::vector<size_t> faceOwnerCells;
    faceOwnerCells.reserve(allFaces.size());
    
    // Add all faces as polygons
    for (const auto& face : allFaces) 
    {
        // Create a polygon for this face
        vtkSmartPointer<vtkPolygon> polygon = 
            vtkSmartPointer<vtkPolygon>::New();

        polygon->GetPointIds()->SetNumberOfIds(face.nodeIndices().size());
        
        // Set the node indices for this polygon
        for (size_t i = 0; i < face.nodeIndices().size(); ++i)
        {
            polygon->GetPointIds()->SetId(i, face.nodeIndices()[i]);
        }
        
        // Add the polygon to the cell array
        polygons->InsertNextCell(polygon);
        
        // Store the owner cell index for this face
        faceOwnerCells.push_back(face.ownerCell());
    }
    
    // Set the polygons to the polyData
    polyData->SetPolys(polygons);
    
    // Add scalar cell fields to the polyData
    // Since PolyData uses faces as cells, we need to map cell data to face data
    for (const auto& [fieldName, scalarField] : scalarCellFields)
    {
        if (!scalarField) continue;

        // Create a vtkDoubleArray for this field
        vtkSmartPointer<vtkDoubleArray> dataArray =
            vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(allFaces.size());

        // Map cell data to face data using the owner cell index
        for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
        {
            size_t ownerCellIdx = faceOwnerCells[faceIdx];

            // Check if the owner cell index is valid
            if (ownerCellIdx < scalarField->size())
            {
                dataArray->SetValue(faceIdx, (*scalarField)[ownerCellIdx]);
            }
            else
            {
                // Set a default value if the owner cell index is out of bounds
                dataArray->SetValue(faceIdx, 0.0);
                std::cerr   << "Warning: Owner cell index " << ownerCellIdx
                            << " is out of bounds for field " << fieldName
                            << std::endl;
            }
        }

        // Add the data array to the polyData's cell data
        polyData->GetCellData()->AddArray(dataArray);
    }

    // Add vector cell fields to the polyData
    for (const auto& [fieldName, vectorField] : vectorCellFields)
    {
        if (!vectorField) continue;

        // Create a vtkDoubleArray for this vector field
        vtkSmartPointer<vtkDoubleArray> dataArray =
            vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(3);
        dataArray->SetNumberOfTuples(allFaces.size());

        // Map cell data to face data using the owner cell index
        for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
        {
            size_t ownerCellIdx = faceOwnerCells[faceIdx];

            // Check if the owner cell index is valid
            if (ownerCellIdx < vectorField->size())
            {
                const Vector& vec = (*vectorField)[ownerCellIdx];
                double vecData[3] = {vec.x(), vec.y(), vec.z()};
                dataArray->SetTuple(faceIdx, vecData);
            }
            else
            {
                // Set a default value if the owner cell index is out of bounds
                double zeroVec[3] = {0.0, 0.0, 0.0};
                dataArray->SetTuple(faceIdx, zeroVec);
                std::cerr   << "Warning: Owner cell index " << ownerCellIdx
                            << " is out of bounds for field " << fieldName
                            << std::endl;
            }
        }

        // Add the data array to the polyData's cell data
        polyData->GetCellData()->AddArray(dataArray);
    }
    
    // Write the polyData to file using vtkXMLPolyDataWriter
    // Note: The output file should use .vtp extension for VTK PolyData XML format
    // If .vtk extension is used, the file will still be written but may cause confusion
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polyData);
    
    // Set the data mode to ASCII for better readability (optional)
    // writer->SetDataModeToAscii();
    
    // Set the data mode to Binary for better performance and smaller file size
    writer->SetDataModeToBinary();
    
    // Enable compression for smaller file size
    writer->SetCompressorTypeToZLib();
    
    // Write the file
    writer->Write();
    
    std::cout   << "VTK PolyData file written successfully: " 
                << filename << std::endl;

    std::cout   << "  - Number of points: " << allNodes.size() << std::endl;

    std::cout   << "  - Number of polygons: " << allFaces.size() << std::endl;

    std::cout   << "  - Number of scalar fields: " << scalarCellFields.size()
                << std::endl;

    std::cout   << "  - Number of vector fields: " << vectorCellFields.size()
                << std::endl;
}

// Helper function to order hexahedron nodes according to VTK convention
// VTK Hexahedron: nodes 0-3 form bottom quad, nodes 4-7 form top quad
// Uses face topology for robust ordering without axis-alignment assumptions
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
        centroid = centroid + allNodes[nodeId];
    }
    centroid = centroid / 8.0;

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
    auto computeQuadNormal = [&allNodes](const std::vector<size_t>& quadNodes) -> Vector
    {
        Vector v0 = allNodes[quadNodes[0]];
        Vector v1 = allNodes[quadNodes[1]];
        Vector v2 = allNodes[quadNodes[2]];
        Vector edge1 = v1 - v0;
        Vector edge2 = v2 - v0;
        return cross(edge1, edge2);
    };

    auto computeFaceCentroid = [&allNodes](const std::vector<size_t>& faceNodes)
    {
        Vector fc(0.0, 0.0, 0.0);
        for (size_t nid : faceNodes)
        {
            fc = fc + allNodes[nid];
        }
        return fc / static_cast<Scalar>(faceNodes.size());
    };

    Vector bottomNormal = computeQuadNormal(bottomFace);
    Vector topNormal = computeQuadNormal(topFace);
    Vector bottomCentroid = computeFaceCentroid(bottomFace);
    Vector topCentroid = computeFaceCentroid(topFace);

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
    Vector testNormal = computeQuadNormal(orderedBottomNodes);
    Vector outwardDir = computeFaceCentroid(orderedBottomNodes) - centroid;
    if (dot(testNormal, outwardDir) < 0)
    {
        std::reverse(orderedBottomNodes.begin(), orderedBottomNodes.end());
    }

    // Build edge connectivity from side faces to match top to bottom
    // Each side face connects an edge of bottom to an edge of top
    std::map<size_t, std::set<size_t>> nodeConnections;

    for (size_t i = 0; i < faceNodeLists.size(); ++i)
    {
        if (i == bottomFaceIdx || i == topFaceIdx) continue; // Skip top/bottom

        const auto& sideFace = faceNodeLists[i];
        std::vector<size_t> sideBottomNodes, sideTopNodes;

        for (size_t nid : sideFace)
        {
            bool inBottom = std::find(bottomFace.begin(), bottomFace.end(), nid) != bottomFace.end();
            bool inTop = std::find(topFace.begin(), topFace.end(), nid) != topFace.end();

            if (inBottom) sideBottomNodes.push_back(nid);
            if (inTop) sideTopNodes.push_back(nid);
        }

        // Record connections
        for (size_t bn : sideBottomNodes)
        {
            for (size_t tn : sideTopNodes)
            {
                nodeConnections[bn].insert(tn);
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
            const auto& connectedNodes = nodeConnections[bottomNode];

            for (size_t candidate : connectedNodes)
            {
                if (std::find(topFace.begin(), topFace.end(), candidate) != topFace.end())
                {
                    orderedTopNodes[i] = candidate;
                    break;
                }
            }
        }
    }

    // Build final VTK hex ordering
    orderedNodes.reserve(8);
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

// Helper function to order wedge (prism) nodes according to VTK convention
// VTK Wedge: nodes 0,1,2 form bottom triangle, nodes 3,4,5 form top triangle
// Uses actual face topology to determine correct node ordering
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
        centroid = centroid + allNodes[nodeId];
    }
    centroid = centroid / 6.0;

    // Compute face centroids and normals for triangular faces
    auto computeFaceCentroid = [&allNodes](const std::vector<size_t>& faceNodes)
    {
        Vector fc(0.0, 0.0, 0.0);
        for (size_t nid : faceNodes)
        {
            fc = fc + allNodes[nid];
        }
        return fc / static_cast<Scalar>(faceNodes.size());
    };

    auto computeTriangleNormal = [&allNodes](const std::vector<size_t>& triNodes) -> Vector
    {
        Vector v0 = allNodes[triNodes[0]];
        Vector v1 = allNodes[triNodes[1]];
        Vector v2 = allNodes[triNodes[2]];
        Vector edge1 = v1 - v0;
        Vector edge2 = v2 - v0;
        return cross(edge1, edge2); // Not normalized
    };

    Vector tri0Centroid = computeFaceCentroid(triangularFaces[0]);
    Vector tri1Centroid = computeFaceCentroid(triangularFaces[1]);

    Vector tri0Normal = computeTriangleNormal(triangularFaces[0]);
    Vector tri1Normal = computeTriangleNormal(triangularFaces[1]);

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
    // Two nodes are from bottom triangle (forming an edge), two from top triangle (forming an edge)
    // The edges are connected by the quad face
    std::map<size_t, size_t> bottomToTopMap;

    for (const auto& quadNodes : quadFaces)
    {
        if (quadNodes.size() != 4) continue;

        // Classify nodes as bottom or top
        std::vector<size_t> quadBottomNodes, quadTopNodes;

        for (size_t nid : quadNodes)
        {
            bool inBottom = std::find(bottomTri.begin(), bottomTri.end(), nid) != bottomTri.end();
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
        // In the quad [n0, n1, n2, n3], edges are: n0-n1, n1-n2, n2-n3, n3-n0

        // Check all 4 edges of the quad to find which bottom and top nodes are connected
        for (size_t i = 0; i < 4; ++i)
        {
            size_t n0 = quadNodes[i];
            size_t n1 = quadNodes[(i + 1) % 4];

            bool n0_bottom = std::find(bottomTri.begin(), bottomTri.end(), n0) != bottomTri.end();
            bool n1_bottom = std::find(bottomTri.begin(), bottomTri.end(), n1) != bottomTri.end();

            // If this edge connects a bottom node to a top node, record the connection
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

    // Order bottom triangle nodes consistently (counter-clockwise when viewed from outside)
    std::vector<size_t> orderedBottomNodes = bottomTri;

    // Ensure proper winding by checking normal direction
    Vector testNormal = computeTriangleNormal(orderedBottomNodes);
    Vector outwardDir = computeFaceCentroid(orderedBottomNodes) - centroid;
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
    std::set<size_t> mappedTopNodes(orderedTopNodes.begin(), orderedTopNodes.end());
    std::set<size_t> expectedTopNodes(topTri.begin(), topTri.end());

    if (mappedTopNodes != expectedTopNodes)
    {
        // Mapping is incomplete or incorrect - return empty to trigger fallback
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

// Helper function to order pyramid nodes according to VTK convention
// VTK Pyramid: nodes 0,1,2,3 form quad base, node 4 is apex
std::vector<vtkIdType> orderPyramidNodes
(
    const std::vector<std::vector<size_t>>& faceNodeLists,
    const std::set<size_t>& uniqueNodes,
    const std::vector<Vector>& allNodes
)
{
    std::vector<vtkIdType> orderedNodes;

    if (uniqueNodes.size() != 5)
    {
        return orderedNodes;
    }

    std::vector<size_t> nodeList(uniqueNodes.begin(), uniqueNodes.end());

    // Find apex (node with max z-coordinate, assuming pyramid points up)
    size_t apexNode = nodeList[0];
    Scalar maxZ = allNodes[nodeList[0]].z();
    for (size_t i = 1; i < 5; ++i)
    {
        if (allNodes[nodeList[i]].z() > maxZ)
        {
            maxZ = allNodes[nodeList[i]].z();
            apexNode = nodeList[i];
        }
    }

    // Collect base nodes (all except apex)
    std::vector<size_t> baseNodes;
    for (size_t nodeId : nodeList)
    {
        if (nodeId != apexNode)
        {
            baseNodes.push_back(nodeId);
        }
    }

    if (baseNodes.size() != 4)
    {
        return orderedNodes;
    }

    // Order base nodes by angle around base centroid
    Vector baseCentroid(0.0, 0.0, 0.0);
    for (size_t nodeId : baseNodes)
    {
        baseCentroid = baseCentroid + allNodes[nodeId];
    }
    baseCentroid = baseCentroid / 4.0;

    std::sort(baseNodes.begin(), baseNodes.end(),
        [&allNodes, &baseCentroid](size_t a, size_t b)
        {
            Vector va = allNodes[a] - baseCentroid;
            Vector vb = allNodes[b] - baseCentroid;
            return std::atan2(va.y(), va.x()) < std::atan2(vb.y(), vb.x());
        });

    // Build ordered list: 4 base nodes + apex
    orderedNodes.reserve(5);
    for (size_t nodeId : baseNodes)
    {
        orderedNodes.push_back(static_cast<vtkIdType>(nodeId));
    }
    orderedNodes.push_back(static_cast<vtkIdType>(apexNode));

    return orderedNodes;
}

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
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

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
    size_t numTets = 0, numHexes = 0, numWedges = 0, numPyramids = 0, numConvex = 0;
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
            // Tetrahedron - VTK ordering: any vertex, then 3 neighbors forming base
            cellType = VTK_TETRA;
            numTets++;
            orderedNodes.reserve(4);
            for (size_t nodeId : uniqueNodeSet)
            {
                orderedNodes.push_back(static_cast<vtkIdType>(nodeId));
            }
        }
        else if (numNodes == 8 && numFaces == 6)
        {
            // Hexahedron - need to find proper ordering
            cellType = VTK_HEXAHEDRON;
            numHexes++;
            orderedNodes = orderHexahedronNodes(faceNodeLists, uniqueNodeSet, allNodes);
        }
        else if (numNodes == 6 && numFaces == 5)
        {
            // Wedge (prism)
            cellType = VTK_WEDGE;
            numWedges++;
            orderedNodes = orderWedgeNodes(faceNodeLists, uniqueNodeSet, allNodes);

            if (orderedNodes.empty())
            {
                numWedgeOrderingFailures++;
                // Fallback to convex point set
                cellType = VTK_CONVEX_POINT_SET;
                for (size_t nodeId : uniqueNodeSet)
                {
                    orderedNodes.push_back(static_cast<vtkIdType>(nodeId));
                }
            }
        }
        else if (numNodes == 5 && numFaces == 5)
        {
            // Pyramid
            cellType = VTK_PYRAMID;
            numPyramids++;
            orderedNodes = orderPyramidNodes(faceNodeLists, uniqueNodeSet, allNodes);
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
            vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
            for (size_t nodeId : uniqueNodeSet)
            {
                cellPointIds->InsertNextId(static_cast<vtkIdType>(nodeId));
            }
            unstructuredGrid->InsertNextCell(VTK_CONVEX_POINT_SET, cellPointIds);
        }
        else
        {
            unstructuredGrid->InsertNextCell(cellType,
                                            static_cast<vtkIdType>(orderedNodes.size()),
                                            orderedNodes.data());
        }
    }

    // Add scalar cell fields
    for (const auto& [fieldName, scalarField] : scalarCellFields)
    {
        if (!scalarField) continue;

        vtkSmartPointer<vtkDoubleArray> dataArray =
            vtkSmartPointer<vtkDoubleArray>::New();

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

    std::cout   << "VTK UnstructuredGrid file written successfully: "
                << filename << std::endl;

    std::cout   << "  - Number of points: " << allNodes.size() << std::endl;

    std::cout   << "  - Number of cells: " << allCells.size() << std::endl;

    std::cout   << "  - Cell types:" << std::endl;
    if (numTets > 0)
        std::cout << "    - Tetrahedra: " << numTets << std::endl;
    if (numHexes > 0)
        std::cout << "    - Hexahedra: " << numHexes << std::endl;
    if (numWedges > 0)
        std::cout << "    - Wedges (prisms): " << numWedges << std::endl;
    if (numPyramids > 0)
        std::cout << "    - Pyramids: " << numPyramids << std::endl;
    if (numConvex > 0)
        std::cout << "    - Convex point sets: " << numConvex << std::endl;

    if (numWedgeOrderingFailures > 0)
    {
        std::cout   << "  - WARNING: " << numWedgeOrderingFailures
                    << " wedge cells failed proper ordering (using convex fallback)"
                    << std::endl;
    }

    std::cout   << "  - Number of scalar fields: " << scalarCellFields.size()
                << std::endl;

    std::cout   << "  - Number of vector fields: " << vectorCellFields.size()
                << std::endl;
}

void writePVDTimeSeriesHeader
(
    const std::string& pvdFilename
)
{
    std::ofstream pvdFile(pvdFilename);

    if (!pvdFile.is_open())
    {
        throw std::runtime_error("Failed to open PVD file: " + pvdFilename);
    }

    pvdFile << "<?xml version=\"1.0\"?>\n";
    pvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvdFile << "  <Collection>\n";
    pvdFile << "  </Collection>\n";
    pvdFile << "</VTKFile>\n";

    pvdFile.close();

    std::cout << "PVD time series header written: " << pvdFilename << std::endl;
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
        throw std::runtime_error("Failed to open PVD file for reading: " + pvdFilename);
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
        throw std::runtime_error("Failed to open PVD file for writing: " + pvdFilename);
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

    std::cout << "Added timestep " << timeValue << " to PVD file" << std::endl;
}

} // namespace VtkWriter
