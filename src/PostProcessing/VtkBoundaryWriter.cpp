/******************************************************************************
 * @file VtkBoundaryWriter.cpp
 * @brief Implementation of wall boundary PolyData writer
 *****************************************************************************/

// Own header
#include "VtkBoundaryWriter.h"

// VTK includes
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

// STL includes
#include <iostream>
#include <unordered_map>
#include <vector>

// Header includes
#include "BoundaryPatch.h"


namespace VTK
{

void writeWallBoundaryData
(
    const std::string& filename,
    const Mesh& mesh,
    const std::map<std::string,
    const FaceData<Scalar>*>& scalarFaceFields,
    bool debug
)
{
    const std::span<const Vector> allNodes = mesh.nodes();
    const std::span<const Face> allFaces = mesh.faces();

    // Collect wall boundary faces
    std::vector<size_t> wallFaceIndices;

    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
    {
        const auto& face = allFaces[faceIdx];
        if (!face.isBoundary()) continue;

        const auto& patch = face.patch();

        if (patch.has_value() && patch->get().type() == PatchType::WALL)
        {
            wallFaceIndices.push_back(faceIdx);
        }
    }

    if (wallFaceIndices.empty())
    {
        std::cout
            << "No wall boundary faces found. "
            << "Skipping wall export." << '\n';
        return;
    }

    // Build node remapping (global -> local)
    std::unordered_map<size_t, vtkIdType> nodeMap;
    vtkIdType localId = 0;

    for (size_t faceIdx : wallFaceIndices)
    {
        for (size_t nodeIdx : allFaces[faceIdx].nodeIndices())
        {
            if (nodeMap.find(nodeIdx) == nodeMap.end())
            {
                nodeMap[nodeIdx] = localId++;
            }
        }
    }

    // Create VTK PolyData
    vtkSmartPointer<vtkPolyData> polyData =
        vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    points->SetNumberOfPoints(static_cast<vtkIdType>(nodeMap.size()));

    for (const auto& [globalIdx, localIdx] : nodeMap)
    {
        const Vector& node = allNodes[globalIdx];
        points->SetPoint
        (
            localIdx,
            static_cast<double>(node.x()),
            static_cast<double>(node.y()),
            static_cast<double>(node.z())
        );
    }

    polyData->SetPoints(points);

    // Add wall face polygons
    vtkSmartPointer<vtkCellArray> cells =
        vtkSmartPointer<vtkCellArray>::New();

    for (size_t faceIdx : wallFaceIndices)
    {
        const auto nodeIndices = allFaces[faceIdx].nodeIndices();
        const vtkIdType npts =
            static_cast<vtkIdType>(nodeIndices.size());

        std::vector<vtkIdType> pts(static_cast<size_t>(npts));

        for (vtkIdType j = 0; j < npts; ++j)
        {
            pts[static_cast<size_t>(j)] = nodeMap[nodeIndices[static_cast<size_t>(j)]];
        }

        cells->InsertNextCell(npts, pts.data());
    }

    polyData->SetPolys(cells);

    // Attach face-centered scalar fields
    for (const auto& [fieldName, faceField] : scalarFaceFields)
    {
        if (!faceField) continue;

        vtkSmartPointer<vtkDoubleArray> dataArray =
            vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(static_cast<vtkIdType>(wallFaceIndices.size()));

        for (size_t i = 0; i < wallFaceIndices.size(); ++i)
        {
            const size_t faceIdx = wallFaceIndices[i];
            const double value = static_cast<double>((*faceField)[faceIdx]);
            dataArray->SetValue(static_cast<vtkIdType>(i), value);
        }

        polyData->GetCellData()->AddArray(dataArray);
    }

    // Write .vtp file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    writer->SetFileName(filename.c_str());
    writer->SetInputData(polyData);
    writer->SetDataModeToBinary();
    writer->SetCompressorTypeToZLib();
    writer->Write();

    std::cout
        << "VTK PolyData wall boundary file written: "
        << filename << '\n';

    if (debug)
    {
        std::cout
            << "  - Wall faces: "
            << wallFaceIndices.size() << '\n';

        std::cout
            << "  - Wall nodes: "
            << nodeMap.size() << '\n';

        std::cout
            << "  - Scalar fields: "
            << scalarFaceFields.size() << '\n';
    }
}

} // namespace VTK
