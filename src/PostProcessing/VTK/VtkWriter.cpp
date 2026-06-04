/******************************************************************************
 * @file VtkWriter.cpp
 * @brief Implementation of VTK UnstructuredGrid writer
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "VtkWriter.h"

// External library headers
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCellValidator.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

// Standard library headers
#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// Project headers
#include "ErrorHandler.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

namespace
{

Vector newellNormal
(
    IndexListRef faceNodes,
    NodeListRef allNodes
)
{
    Vector normal;

    for (Index i = 0; i < faceNodes.size(); ++i)
    {
        const Vector& current = allNodes[faceNodes[i]];
        const Vector& next = allNodes[faceNodes[(i + 1) % faceNodes.size()]];

        normal += Vector
        (
            (current.y() - next.y()) * (current.z() + next.z()),
            (current.z() - next.z()) * (current.x() + next.x()),
            (current.x() - next.x()) * (current.y() + next.y())
        );
    }

    return normal;
}


void appendUniqueNodeId
(
    std::vector<vtkIdType>& uniquePointIds,
    Index nodeIdx
)
{
    const vtkIdType vtkNodeIdx = static_cast<vtkIdType>(nodeIdx);

    if
    (
        std::find
        (
            uniquePointIds.begin(),
            uniquePointIds.end(),
            vtkNodeIdx
        ) == uniquePointIds.end()
    )
    {
        uniquePointIds.push_back(vtkNodeIdx);
    }
}


void appendValidationStateName
(
    Message& result,
    const vtkCellValidator::State state,
    const vtkCellValidator::State bit,
    const char* const name
)
{
    if ((state & bit) == bit)
    {
        if (!result.empty())
        {
            result += '|';
        }

        result += name;
    }
}


Message cellValidationStateName(vtkCellValidator::State state)
{
    if (state == vtkCellValidator::Valid)
    {
        return "Valid";
    }

    Message result;

    appendValidationStateName
    (
        result,
        state,
        vtkCellValidator::WrongNumberOfPoints,
        "WrongNumberOfPoints"
    );
    appendValidationStateName
    (
        result,
        state,
        vtkCellValidator::IntersectingEdges,
        "IntersectingEdges"
    );
    appendValidationStateName
    (
        result,
        state,
        vtkCellValidator::IntersectingFaces,
        "IntersectingFaces"
    );
    appendValidationStateName
    (
        result,
        state,
        vtkCellValidator::NoncontiguousEdges,
        "NoncontiguousEdges"
    );
    appendValidationStateName
    (
        result,
        state,
        vtkCellValidator::Nonconvex,
        "Nonconvex"
    );
    appendValidationStateName
    (
        result,
        state,
        vtkCellValidator::FacesAreOrientedIncorrectly,
        "FacesAreOrientedIncorrectly"
    );

    return result.empty() ? "Unknown" : result;
}


bool hasState
(
    const vtkCellValidator::State state,
    const vtkCellValidator::State bit
)
{
    return (state & bit) == bit;
}


void validateCellsForDebug(vtkUnstructuredGrid* unstructuredGrid)
{
    vtkSmartPointer<vtkGenericCell> genericCell =
        vtkSmartPointer<vtkGenericCell>::New();

    // vtkCellValidator discounts the shared vertices of adjacent edges and
    // faces by comparing parametric coordinates (in [0, 1]) against this
    // tolerance. FLT_EPSILON is the value vtkCellValidator itself defaults
    // to; machine epsilon (smallValue in double precision) is roughly 5e8
    // times tighter, so ordinary round-off in vtkLine::Intersection reads
    // shared vertices as genuine crossings and spuriously flags warped cells
    // as IntersectingEdges / FacesAreOrientedIncorrectly.
    const double tolerance =
        static_cast<double>(std::numeric_limits<float>::epsilon());

    Count invalidCellCount = 0;
    Count nonConvexCellCount = 0;
    std::vector<Message> invalidSamples;
    constexpr Count maxSamples = 8;

    const vtkIdType numCells = unstructuredGrid->GetNumberOfCells();

    for (vtkIdType cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        unstructuredGrid->GetCell(cellIdx, genericCell);

        const vtkCellValidator::State state =
            vtkCellValidator::Check(genericCell.GetPointer(), tolerance);

        if (state == vtkCellValidator::Valid)
        {
            continue;
        }

        // Nonconvex is a mesh-quality property, not a defect in how this
        // writer built the cell: a non-convex polyhedron is a valid VTK
        // cell that renders correctly. Only genuine structural problems
        // count as invalid output, so a cell flagged solely as Nonconvex is
        // reported separately and does not trigger the warning.
        const bool structurallyInvalid =
            hasState(state, vtkCellValidator::WrongNumberOfPoints)
         || hasState(state, vtkCellValidator::IntersectingEdges)
         || hasState(state, vtkCellValidator::IntersectingFaces)
         || hasState(state, vtkCellValidator::NoncontiguousEdges)
         || hasState(state, vtkCellValidator::FacesAreOrientedIncorrectly);

        if (!structurallyInvalid)
        {
            ++nonConvexCellCount;
            continue;
        }

        ++invalidCellCount;

        if (invalidSamples.size() < maxSamples)
        {
            std::ostringstream sample;
            sample
                << cellIdx << '=' << cellValidationStateName(state);
            invalidSamples.push_back(sample.str());
        }
    }

    if (nonConvexCellCount > 0)
    {
        std::cout
            << "VTK cell validation: " << nonConvexCellCount
            << " non-convex cell(s) (mesh quality; valid for output).\n";
    }

    if (invalidCellCount == 0)
    {
        return;
    }

    std::ostringstream message;
    message
        << "VTK cell validation reported " << invalidCellCount
        << " structurally invalid cell(s). Output will still be written.";

    if (!invalidSamples.empty())
    {
        message << " First samples: ";

        for (Index i = 0; i < invalidSamples.size(); ++i)
        {
            if (i > 0)
            {
                message << ", ";
            }

            message << invalidSamples[i];
        }
    }

    Warning(message.str());
}

} // namespace


void writeVtkUnstructuredGrid
(
    const FilePath& filename,
    const Mesh& mesh,
    const ScalarFieldMap& scalarFields,
    const VectorFieldMap& vectorFields,
    bool debug
)
{
    const NodeListRef allNodes = mesh.nodes();
    const CellListRef allCells = mesh.cells();
    const FaceListRef allFaces = mesh.faces();

    // Create vtkUnstructuredGrid
    vtkSmartPointer<vtkUnstructuredGrid>
    unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Create vtkPoints from nodes
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(static_cast<vtkIdType>(allNodes.size()));

    for (Index nodeIdx = 0; nodeIdx < allNodes.size(); ++nodeIdx)
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
    unstructuredGrid->Allocate(static_cast<vtkIdType>(allCells.size()));

    // Add cells as polyhedra using the mesh's face-based topology directly.
    for (Index cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        const Cell& cell = allCells[cellIdx];
        const auto faceIndices = cell.faceIndices();

        std::vector<vtkIdType> uniquePointIds;
        std::vector<vtkIdType> faceStream;
        uniquePointIds.reserve(faceIndices.size() * 4);
        faceStream.reserve(faceIndices.size() * 5);

        for (Index faceIdx : faceIndices)
        {
            const Face& face = allFaces[faceIdx];
            const auto nodeIndices = face.nodeIndices();

            if (nodeIndices.size() < 3)
            {
                FatalError
                (
                    "Cannot write VTK polyhedron cell "
                  + std::to_string(cellIdx) + ": face "
                  + std::to_string(faceIdx)
                  + " has fewer than three nodes."
                );
            }

            for (Index nodeIdx : nodeIndices)
            {
                appendUniqueNodeId(uniquePointIds, nodeIdx);
            }

            std::vector<vtkIdType> orientedFaceNodes;
            orientedFaceNodes.reserve(nodeIndices.size());

            for (Index nodeIdx : nodeIndices)
            {
                orientedFaceNodes.push_back(static_cast<vtkIdType>(nodeIdx));
            }

            const Vector faceNormal = newellNormal(nodeIndices, allNodes);
            const Vector outwardDirection =
                face.centroid() - cell.centroid();

            if (dot(faceNormal, outwardDirection) < S(0.0))
            {
                std::reverse
                (
                    orientedFaceNodes.begin(),
                    orientedFaceNodes.end()
                );
            }

            faceStream.push_back
            (
                static_cast<vtkIdType>(orientedFaceNodes.size())
            );
            faceStream.insert
            (
                faceStream.end(),
                orientedFaceNodes.begin(),
                orientedFaceNodes.end()
            );
        }

        const vtkIdType insertedCell = unstructuredGrid->InsertNextCell
        (
            VTK_POLYHEDRON,
            static_cast<vtkIdType>(uniquePointIds.size()),
            uniquePointIds.data(),
            static_cast<vtkIdType>(faceIndices.size()),
            faceStream.data()
        );

        if (insertedCell < 0)
        {
            FatalError
            (
                "Failed to insert VTK polyhedron cell "
              + std::to_string(cellIdx) + "."
            );
        }
    }

    // Add scalar cell fields
    for (const auto& [fieldName, scalarField] : scalarFields)
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

        for (Index cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            const double value =
                static_cast<double>((*scalarField)[cellIdx]);
            dataArray->SetValue(static_cast<vtkIdType>(cellIdx), value);
        }

        unstructuredGrid->GetCellData()->AddArray(dataArray);
    }

    // Add vector cell fields (supplied as three scalar component fields)
    for (const auto& [fieldName, components] : vectorFields)
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

        for (Index cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
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

    if (debug)
    {
        validateCellsForDebug(unstructuredGrid);
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

    if (debug)
    {
        std::cout
            << "VTK UnstructuredGrid file written: "
            << filename << '\n';

        std::cout
            << "  - Number of points: " << allNodes.size() << '\n';

        std::cout
            << "  - Number of cells: " << allCells.size() << '\n';

        std::cout
            << "  - Cell type: VTK_POLYHEDRON" << '\n';

        std::cout
            << "  - Number of scalar fields: "
            << scalarFields.size() << '\n';

        std::cout
            << "  - Number of vector fields: "
            << vectorFields.size() << '\n';
    }
}

} // namespace VTK
