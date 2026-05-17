/******************************************************************************
 * @file VtkWriter.cpp
 * @brief Implementation of VTK UnstructuredGrid writer
 *****************************************************************************/

#include "VtkWriter.h"

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include "ErrorHandler.h"
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
    std::array<const ScalarField*, 3>>& vectorCellFields,
    bool debug
)
{
    const std::span<const Vector> allNodes = mesh.nodes();
    const std::span<const Cell> allCells = mesh.cells();
    const std::span<const Face> allFaces = mesh.faces();

    // Create vtkUnstructuredGrid
    vtkSmartPointer<vtkUnstructuredGrid>
    unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Create vtkPoints from nodes
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(static_cast<vtkIdType>(allNodes.size()));

    for (size_t nodeIdx = 0; nodeIdx < allNodes.size(); ++nodeIdx)
    {
        const Vector& node = allNodes[nodeIdx];
        points->SetPoint
        (
            static_cast<vtkIdType>(nodeIdx),
            static_cast<double>(node.x()),
            static_cast<double>(node.y()),
            static_cast<double>(node.z())
        );
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
        const auto faceIndices = cell.faceIndices();

        // Collect all unique nodes and build connectivity
        std::unordered_set<size_t> uniqueNodeSet;
        std::vector<std::vector<size_t>> faceNodeLists;
        faceNodeLists.reserve(faceIndices.size());

        for (size_t faceIdx : faceIndices)
        {
            const Face& face = allFaces[faceIdx];
            const auto nodeIndices = face.nodeIndices();
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

        if (scalarField->size() != allCells.size())
        {
            FatalError
            (
                "VTK scalar field '" + fieldName + "' has "
              + std::to_string(scalarField->size())
              + " values but the mesh has "
              + std::to_string(allCells.size()) + " cells."
            );
        }

        vtkSmartPointer<vtkDoubleArray>
        dataArray = vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(static_cast<vtkIdType>(allCells.size()));

        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            const double value =
                static_cast<double>((*scalarField)[cellIdx]);
            dataArray->SetValue(static_cast<vtkIdType>(cellIdx), value);
        }

        unstructuredGrid->GetCellData()->AddArray(dataArray);
    }

    // Add vector cell fields (supplied as three scalar component fields)
    for (const auto& [fieldName, components] : vectorCellFields)
    {
        if (!components[0] || !components[1] || !components[2]) continue;

        const ScalarField& cx = *components[0];
        const ScalarField& cy = *components[1];
        const ScalarField& cz = *components[2];

        if
        (
            cx.size() != allCells.size()
         || cy.size() != allCells.size()
         || cz.size() != allCells.size()
        )
        {
            FatalError
            (
                "VTK vector field '" + fieldName
              + "' has component sizes that do not match the mesh cell "
                "count " + std::to_string(allCells.size()) + "."
            );
        }

        vtkSmartPointer<vtkDoubleArray> dataArray =
            vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(3);
        dataArray->SetNumberOfTuples(static_cast<vtkIdType>(allCells.size()));

        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            const double vecData[3] =
            {
                static_cast<double>(cx[cellIdx]),
                static_cast<double>(cy[cellIdx]),
                static_cast<double>(cz[cellIdx])
            };
            dataArray->SetTuple(static_cast<vtkIdType>(cellIdx), vecData);
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

    if (writer->Write() == 0)
    {
        FatalError("Failed to write VTK UnstructuredGrid file: " + filename);
    }

    std::cout
        << "VTK UnstructuredGrid file written: "
        << filename << '\n';

    if (debug)
    {
        std::cout
            << "  - Number of points: " << allNodes.size() << '\n';

        std::cout
            << "  - Number of cells: " << allCells.size() << '\n';

        std::cout
            << "  - Cell types:" << '\n';

        if (numTets > 0)
        {
            std::cout
                << "    - Tetrahedra: " << numTets << '\n';
        }

        if (numHexes > 0)
        {
            std::cout
                << "    - Hexahedra: " << numHexes << '\n';
        }

        if (numWedges > 0)
        {
            std::cout
                << "    - Wedges (prisms): " << numWedges << '\n';
        }

        if (numPyramids > 0)
        {
            std::cout
                << "    - Pyramids: " << numPyramids << '\n';
        }

        if (numConvex > 0)
        {
            std::cout
                << "    - Convex point sets: " << numConvex << '\n';
        }

        if (numWedgeOrderingFailures > 0)
        {
            std::cout
                << "  - WARNING: " << numWedgeOrderingFailures
                << " wedge cells failed proper ordering" << '\n';
        }

        std::cout
            << "  - Number of scalar fields: "
            << scalarCellFields.size() << '\n';

        std::cout
            << "  - Number of vector fields: "
            << vectorCellFields.size() << '\n';
    }
}

} // namespace VTK
