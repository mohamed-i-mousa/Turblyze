/******************************************************************************
 * @file VtkWriter.cpp
 * @brief Implementation of VTK UnstructuredGrid writer
 *****************************************************************************/

// Own header
#include "VtkWriter.h"

// VTK includes
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

// STL includes
#include <iostream>
#include <unordered_set>
#include <vector>

// Header includes
#include "VtkCellOrdering.h"


namespace VTK
{

void writeVtkUnstructuredGrid
(
    const std::string& filename,
    const Mesh& mesh,
    const std::map<std::string,
    const ScalarField*>& scalarCellFields,
    const std::map<std::string,
    const VectorField*>& vectorCellFields,
    bool debug
)
{
    std::span<const Vector> allNodes = mesh.nodes();
    std::span<const Cell> allCells = mesh.cells();
    std::span<const Face> allFaces = mesh.faces();

    // Create vtkUnstructuredGrid
    vtkSmartPointer<vtkUnstructuredGrid>
    unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Create vtkPoints from nodes
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(static_cast<vtkIdType>(allNodes.size()));

    for (size_t nodeIdx = 0; nodeIdx < allNodes.size(); ++nodeIdx)
    {
        points->SetPoint(static_cast<vtkIdType>(nodeIdx), allNodes[nodeIdx].x(), allNodes[nodeIdx].y(), allNodes[nodeIdx].z());
    }

    unstructuredGrid->SetPoints(points);

    // Add cells to unstructured grid with proper topology preservation
    // Build a node connectivity graph to determine proper VTK node ordering

    // Cell type statistics for diagnostics
    size_t numTets = 0;
    size_t numHexes = 0;
    size_t numWedges = 0;
    size_t numPyramids = 0;
    size_t numConvex = 0;
    size_t numWedgeOrderingFailures = 0;

    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        const Cell& cell = allCells[cellIdx];
        const auto& faceIndices = cell.faceIndices();

        // Collect all unique nodes and build connectivity
        std::unordered_set<size_t> uniqueNodeSet;
        std::vector<std::vector<size_t>> faceNodeLists;
        faceNodeLists.reserve(faceIndices.size());

        for (size_t faceIdx : faceIndices)
        {
            const Face& face = allFaces[faceIdx];
            const auto& nodeIndices = face.nodeIndices();
            uniqueNodeSet.insert(nodeIndices.begin(), nodeIndices.end());
            faceNodeLists.emplace_back(nodeIndices.begin(), nodeIndices.end());
        }

        const size_t numNodes = uniqueNodeSet.size();
        const size_t numFaces = faceIndices.size();

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
            for (size_t nodeIdx : uniqueNodeSet)
            {
                orderedNodes.push_back(static_cast<vtkIdType>(nodeIdx));
            }
        }

        // Insert cell with ordered nodes
        if (orderedNodes.empty())
        {
            // Fallback: use arbitrary ordering
            orderedNodes.reserve(numNodes);
            for (size_t nodeIdx : uniqueNodeSet)
            {
                orderedNodes.push_back(static_cast<vtkIdType>(nodeIdx));
            }
            cellType = VTK_CONVEX_POINT_SET;
        }

        unstructuredGrid->InsertNextCell
        (
            cellType,
            static_cast<vtkIdType>(orderedNodes.size()),
            orderedNodes.data()
        );
    }

    // Add scalar cell fields
    for (const auto& [fieldName, scalarField] : scalarCellFields)
    {
        if (!scalarField) continue;

        vtkSmartPointer<vtkDoubleArray>
        dataArray = vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(static_cast<vtkIdType>(allCells.size()));

        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            if (cellIdx < scalarField->size())
            {
                dataArray->SetValue(static_cast<vtkIdType>(cellIdx), (*scalarField)[cellIdx]);
            }
            else
            {
                dataArray->SetValue(static_cast<vtkIdType>(cellIdx), 0.0);
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
        dataArray->SetNumberOfTuples(static_cast<vtkIdType>(allCells.size()));

        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            if (cellIdx < vectorField->size())
            {
                const Vector& vec = (*vectorField)[cellIdx];
                double vecData[3] = {vec.x(), vec.y(), vec.z()};
                dataArray->SetTuple(static_cast<vtkIdType>(cellIdx), vecData);
            }
            else
            {
                double zeroVec[3] = {0.0, 0.0, 0.0};
                dataArray->SetTuple(static_cast<vtkIdType>(cellIdx), zeroVec);
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
        << "VTK UnstructuredGrid file written: "
        << filename << std::endl;

    if (debug)
    {
        std::cout
            << "  - Number of points: "
            << allNodes.size() << std::endl;

        std::cout
            << "  - Number of cells: "
            << allCells.size() << std::endl;

        std::cout
            << "  - Cell types:" << std::endl;

        if (numTets > 0)
        {
            std::cout
                << "    - Tetrahedra: "
                << numTets << std::endl;
        }

        if (numHexes > 0)
        {
            std::cout
                << "    - Hexahedra: "
                << numHexes << std::endl;
        }

        if (numWedges > 0)
        {
            std::cout
                << "    - Wedges (prisms): "
                << numWedges << std::endl;
        }

        if (numPyramids > 0)
        {
            std::cout
                << "    - Pyramids: "
                << numPyramids << std::endl;
        }

        if (numConvex > 0)
        {
            std::cout
                << "    - Convex point sets: "
                << numConvex << std::endl;
        }

        if (numWedgeOrderingFailures > 0)
        {
            std::cout
                << "  - WARNING: "
                << numWedgeOrderingFailures
                << " wedge cells failed proper ordering"
                << std::endl;
        }

        std::cout
            << "  - Number of scalar fields: "
            << scalarCellFields.size()
            << std::endl;

        std::cout
            << "  - Number of vector fields: "
            << vectorCellFields.size()
            << std::endl;
    }
}

} // namespace VTK
