/******************************************************************************
 * @file VtkBoundaryWriter.cpp
 * @brief Implementation of boundary PolyData writer
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "VtkBoundaryWriter.h"

// Standard library headers
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

// External library headers
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLPolyDataWriter.h>

// Project headers
#include "BoundaryPatch.h"
#include "ErrorHandler.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

// ***************************** Internal Helpers *****************************

namespace
{

struct BoundaryFace
{
    size_t faceIdx = 0;
    int patchIdx = 0;
    int patchZoneIdx = 0;
    int patchTypeIdx = 0;
    unsigned char isWall = 0;
};


int checkedInt(size_t value, const std::string& label)
{
    if (value > static_cast<size_t>(std::numeric_limits<int>::max()))
    {
        FatalError(label + " exceeds the range of vtkIntArray.");
    }

    return static_cast<int>(value);
}


vtkIdType localPointIdx
(
    std::unordered_map<size_t, vtkIdType>& nodeMap,
    std::vector<size_t>& localToGlobalNodes,
    const size_t globalNodeIdx
)
{
    const auto it = nodeMap.find(globalNodeIdx);
    if (it != nodeMap.end())
    {
        return it->second;
    }

    const vtkIdType localIdx =
        static_cast<vtkIdType>(localToGlobalNodes.size());

    nodeMap.emplace(globalNodeIdx, localIdx);
    localToGlobalNodes.push_back(globalNodeIdx);
    return localIdx;
}

} // namespace

// **************************** Write Boundary Data ***************************

void writeBoundaryData
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
    const std::span<const BoundaryPatch> allPatches = mesh.patches();

    std::vector<BoundaryFace> boundaryFaces;
    std::unordered_map<size_t, vtkIdType> nodeMap;
    std::vector<size_t> localToGlobalNodes;

    for (size_t patchIdx = 0; patchIdx < allPatches.size(); ++patchIdx)
    {
        const BoundaryPatch& patch = allPatches[patchIdx];
        const int patchIdxInt = checkedInt(patchIdx, "patchIdx");
        const int patchZoneIdx = checkedInt(patch.zoneIdx(), "patchZoneIdx");
        const int patchTypeIdx = static_cast<int>(patch.type());
        const unsigned char isWall =
            static_cast<unsigned char>
            (
                patch.type() == PatchType::wall ? 1 : 0
            );

        for
        (
            size_t faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            if (faceIdx >= allFaces.size())
            {
                FatalError
                (
                    "Boundary patch '" + patch.patchName()
                  + "' references face " + std::to_string(faceIdx)
                  + ", but the mesh has only "
                  + std::to_string(allFaces.size()) + " faces."
                );
            }

            const Face& face = allFaces[faceIdx];

            if (!face.isBoundary())
            {
                FatalError
                (
                    "Boundary patch '" + patch.patchName()
                  + "' references internal face "
                  + std::to_string(faceIdx) + "."
                );
            }

            for (size_t nodeIdx : face.nodeIndices())
            {
                localPointIdx(nodeMap, localToGlobalNodes, nodeIdx);
            }

            boundaryFaces.push_back
            (
                BoundaryFace
                {
                    faceIdx,
                    patchIdxInt,
                    patchZoneIdx,
                    patchTypeIdx,
                    isWall
                }
            );
        }
    }

    if (boundaryFaces.empty())
    {
        std::cout
            << "No boundary faces found. "
            << "Skipping boundary export." << '\n';
        return;
    }

    for (const auto& [fieldName, faceField] : scalarFaceFields)
    {
        if (!faceField) continue;

        if (faceField->size() != allFaces.size())
        {
            FatalError
            (
                "VTK boundary scalar field '" + fieldName + "' has "
              + std::to_string(faceField->size())
              + " values but the mesh has "
              + std::to_string(allFaces.size()) + " faces."
            );
        }
    }

    // Create VTK PolyData
    vtkSmartPointer<vtkPolyData> polyData =
        vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    points->SetDataTypeToDouble();
    points->SetNumberOfPoints
    (
        static_cast<vtkIdType>(localToGlobalNodes.size())
    );

    for (size_t localIdx = 0; localIdx < localToGlobalNodes.size(); ++localIdx)
    {
        const size_t globalIdx = localToGlobalNodes[localIdx];
        const Vector& node = allNodes[globalIdx];
        points->SetPoint
        (
            static_cast<vtkIdType>(localIdx),
            static_cast<double>(node.x()),
            static_cast<double>(node.y()),
            static_cast<double>(node.z())
        );
    }

    polyData->SetPoints(points);

    // Add boundary face polygons
    vtkSmartPointer<vtkCellArray> cells =
        vtkSmartPointer<vtkCellArray>::New();

    for (const BoundaryFace& record : boundaryFaces)
    {
        const auto nodeIndices = allFaces[record.faceIdx].nodeIndices();
        const vtkIdType npts =
            static_cast<vtkIdType>(nodeIndices.size());

        std::vector<vtkIdType> pts(static_cast<size_t>(npts));

        for (vtkIdType j = 0; j < npts; ++j)
        {
            const size_t nodeIdx = nodeIndices[static_cast<size_t>(j)];
            pts[static_cast<size_t>(j)] = nodeMap[nodeIdx];
        }

        cells->InsertNextCell(npts, pts.data());
    }

    polyData->SetPolys(cells);

    vtkSmartPointer<vtkIntArray> patchIdxArray =
        vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> patchZoneIdxArray =
        vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> patchTypeIdxArray =
        vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkUnsignedCharArray> isWallArray =
        vtkSmartPointer<vtkUnsignedCharArray>::New();

    patchIdxArray->SetName("patchIdx");
    patchZoneIdxArray->SetName("patchZoneIdx");
    patchTypeIdxArray->SetName("patchTypeIdx");
    isWallArray->SetName("isWall");

    patchIdxArray->SetNumberOfComponents(1);
    patchZoneIdxArray->SetNumberOfComponents(1);
    patchTypeIdxArray->SetNumberOfComponents(1);
    isWallArray->SetNumberOfComponents(1);

    const vtkIdType numBoundaryFaces =
        static_cast<vtkIdType>(boundaryFaces.size());

    patchIdxArray->SetNumberOfTuples(numBoundaryFaces);
    patchZoneIdxArray->SetNumberOfTuples(numBoundaryFaces);
    patchTypeIdxArray->SetNumberOfTuples(numBoundaryFaces);
    isWallArray->SetNumberOfTuples(numBoundaryFaces);

    for (size_t i = 0; i < boundaryFaces.size(); ++i)
    {
        const BoundaryFace& record = boundaryFaces[i];
        const vtkIdType tupleIdx = static_cast<vtkIdType>(i);

        patchIdxArray->SetValue(tupleIdx, record.patchIdx);
        patchZoneIdxArray->SetValue(tupleIdx, record.patchZoneIdx);
        patchTypeIdxArray->SetValue(tupleIdx, record.patchTypeIdx);
        isWallArray->SetValue(tupleIdx, record.isWall);
    }

    polyData->GetCellData()->AddArray(patchIdxArray);
    polyData->GetCellData()->AddArray(patchZoneIdxArray);
    polyData->GetCellData()->AddArray(patchTypeIdxArray);
    polyData->GetCellData()->AddArray(isWallArray);

    // Attach face-centered scalar fields
    for (const auto& [fieldName, faceField] : scalarFaceFields)
    {
        if (!faceField) continue;

        vtkSmartPointer<vtkDoubleArray> dataArray =
            vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(numBoundaryFaces);

        for (size_t i = 0; i < boundaryFaces.size(); ++i)
        {
            const size_t faceIdx = boundaryFaces[i].faceIdx;
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

    if (writer->Write() == 0)
    {
        FatalError("Failed to write VTK PolyData boundary file: " + filename);
    }

    std::cout
        << "VTK PolyData boundary file written: "
        << filename << '\n';

    if (debug)
    {
        std::cout
            << "  - Boundary faces: "
            << boundaryFaces.size() << '\n';

        std::cout
            << "  - Boundary nodes: "
            << localToGlobalNodes.size() << '\n';

        std::cout
            << "  - Boundary patches: "
            << allPatches.size() << '\n';

        std::cout
            << "  - Scalar fields: "
            << scalarFaceFields.size() << '\n';
    }
}

} // namespace VTK
