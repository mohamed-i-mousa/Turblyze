/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HDF5BoundaryData.cpp
 * @brief Implementation of the VTKHDF boundary data
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "HDF5BoundaryData.h"

// Standard library headers
#include <algorithm>
#include <format>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

// Project headers
#include "BoundaryConditions.h"
#include "BoundaryPatch.h"
#include "Comm.h"
#include "ErrorHandler.h"
#include "HDF5Wrapper.h"
#include "Reduce.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

// ***************************** Internal Helpers *****************************

namespace
{

// Check that a Count value fits in a 32-bit int, for VTKHDF array sizes
int checkedInt(Count value, const Message& label)
{
    if (value > static_cast<Count>(std::numeric_limits<int>::max()))
    {
        FatalError(label + " exceeds the range of a 32-bit VTKHDF array.");
    }

    return static_cast<int>(value);
}


// Write one PolyData topology group with write-once datasets
void writeTopologyGroup
(
    hid_t location,
    const char* name,
    const std::vector<long long>& connectivity,
    const std::vector<long long>& offsets,
    long long numCells,
    const HDF5::SlabLayout& connectivityRows,
    const HDF5::SlabLayout& offsetsRows
)
{
    hid_t group = HDF5::createGroup(location, name);

    const long long numConnectivityIds =
        static_cast<long long>(connectivity.size());

    const HDF5::SlabLayout partRow = HDF5::oneRowPerRank();

    HDF5::writeDataset
    (
        group,
        "NumberOfCells",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numCells,
        partRow,
        0
    );
    HDF5::writeDataset
    (
        group,
        "NumberOfConnectivityIds",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numConnectivityIds,
        partRow,
        0
    );
    HDF5::writeDataset
    (
        group,
        "Connectivity",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        connectivity.data(),
        connectivityRows,
        0
    );
    HDF5::writeDataset
    (
        group,
        "Offsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        offsets.data(),
        offsetsRows,
        0
    );

    HDF5::closeGroup(group);
}

} // namespace

// ************************* Special Member Functions *************************

HDF5BoundaryData::HDF5BoundaryData
(
    FilePath fileName,
    const Mesh& mesh,
    const BoundaryConditions* bc,
    bool debug
)
:
    mesh_{mesh},
    bc_{bc},
    fileName_{std::move(fileName)},
    debug_{debug},
    file_{H5I_INVALID_HID},
    vtkhdfGroup_{H5I_INVALID_HID},
    cellDataGroup_{H5I_INVALID_HID},
    stepsGroup_{H5I_INVALID_HID},
    cellDataOffsetsGroup_{H5I_INVALID_HID},
    geometryWritten_{false},
    finalized_{false},
    stepCount_{0},
    globalNumBoundaryFaces_{0},
    boundaryFaceRowOffset_{0}
{
    file_ = HDF5::createFile(fileName_);
    vtkhdfGroup_ = HDF5::createGroup(file_, "VTKHDF");

    const long long version[2] = {2, 4};
    HDF5::writeIntArrayAttribute(vtkhdfGroup_, "Version", version, 2);
    HDF5::writeStringAttribute(vtkhdfGroup_, "Type", "PolyData");

    cellDataGroup_ = HDF5::createGroup(vtkhdfGroup_, "CellData");
    stepsGroup_ = HDF5::createGroup(vtkhdfGroup_, "Steps");
    cellDataOffsetsGroup_ = HDF5::createGroup(stepsGroup_, "CellDataOffsets");

    HDF5::writeIntScalarAttribute(stepsGroup_, "NSteps", 0);

    // Keep the file closed between writes so it stays readable mid-run
    closeHandles();
}


HDF5BoundaryData::~HDF5BoundaryData() noexcept
{
    finalize();
}

// ****************************** Public Methods ******************************

void HDF5BoundaryData::writeGeometry()
{
    if (geometryWritten_)
    {
        FatalError
        (
            "VTKHDF boundary geometry has already been written: " + fileName_
        );
    }

    openHandles();

    const NodeList& allNodes = mesh_.nodes();
    const FaceList& allFaces = mesh_.faces();
    const PatchList& allPatches = mesh_.patches();

    // Collect boundary patch faces over a compressed local point set
    std::unordered_map<Index, long long> nodeMap;
    IndexList localToGlobalNodes;

    for (Index patchIdx = 0; patchIdx < allPatches.size(); ++patchIdx)
    {
        const BoundaryPatch& patch = allPatches[patchIdx];

        // Inter-rank cuts are not physical boundary
        if (patch.type() == PatchType::processor)
        {
            continue;
        }

        const int patchIdxInt = checkedInt(patchIdx, "patchIdx");
        const int patchZoneIdx = checkedInt(patch.zoneIdx(), "patchZoneIdx");
        const int patchTypeIdx = static_cast<int>(patch.type());
        const unsigned char isWall =
            static_cast<unsigned char>
            (
                (
                    bc_ != nullptr 
                 && bc_->boundaryType(patch, Field::Ux).isWall()
                ) ? 1 : 0
            );

        for
        (
            Index faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            const Face& face = allFaces[faceIdx];

            if (!face.isBoundary())
            {
                FatalError
                (
                    "Boundary patch '" + patch.name()
                  + "' references internal face "
                  + std::to_string(faceIdx) + "."
                );
            }

            for (Index nodeIdx : face.nodeIndices())
            {
                if (nodeMap.find(nodeIdx) == nodeMap.end())
                {
                    nodeMap.emplace
                    (
                        nodeIdx,
                        static_cast<long long>(localToGlobalNodes.size())
                    );
                    localToGlobalNodes.push_back(nodeIdx);
                }
            }

            boundaryFaceIndices_.push_back(faceIdx);
            patchIdx_.push_back(patchIdxInt);
            patchZoneIdx_.push_back(patchZoneIdx);
            patchTypeIdx_.push_back(patchTypeIdx);
            isWall_.push_back(isWall);
        }
    }

    // Only a globally empty boundary is worth a warning
    const Count globalBoundaryFaces = globalSum(boundaryFaceIndices_.size());

    if (globalBoundaryFaces == 0 && Comm::master())
    {
        Warning
        (
            "No boundary faces found; the boundary VTKHDF file will hold "
            "empty geometry: " + fileName_
        );
    }

    // Local point coordinates, written once as double precision
    std::vector<double> points;
    points.reserve(localToGlobalNodes.size() * 3);

    for (Index globalIdx : localToGlobalNodes)
    {
        const Vector& node = allNodes[globalIdx];
        points.push_back(static_cast<double>(node.x()));
        points.push_back(static_cast<double>(node.y()));
        points.push_back(static_cast<double>(node.z()));
    }

    // Polygon connectivity over local point ids
    std::vector<long long> connectivity;
    std::vector<long long> offsets;
    offsets.reserve(boundaryFaceIndices_.size() + 1);
    offsets.push_back(0);

    for (Index faceIdx : boundaryFaceIndices_)
    {
        for (Index nodeIdx : allFaces[faceIdx].nodeIndices())
        {
            connectivity.push_back(nodeMap.at(nodeIdx));
        }

        offsets.push_back(static_cast<long long>(connectivity.size()));
    }

    // Every rank builds these in the same order, empty or not
    const HDF5::SlabLayout pointRows =
        HDF5::distributedRows(localToGlobalNodes.size());
    const HDF5::SlabLayout polygonRows =
        HDF5::distributedRows(boundaryFaceIndices_.size());
    const HDF5::SlabLayout connectivityRows =
        HDF5::distributedRows(connectivity.size());

    const HDF5::SlabLayout partRow = HDF5::oneRowPerRank();

    // Empty on EVERY rank: no data rows, but still one offsets row
    const HDF5::SlabLayout emptyRows = {0, 0, 0};

    const long long numPoints =
        static_cast<long long>(localToGlobalNodes.size());

    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "NumberOfPoints",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numPoints,
        partRow,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "Points",
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        points.data(),
        pointRows,
        3
    );

    // The four PolyData topologies, in VTKHDF order; only Polygons is filled
    const std::vector<long long> emptyConnectivity;
    const std::vector<long long> emptyOffsets(1, 0);

    writeTopologyGroup
    (
        vtkhdfGroup_,
        "Vertices",
        emptyConnectivity,
        emptyOffsets,
        0,
        emptyRows,
        HDF5::perPieceOffsetsRows(emptyRows)
    );
    writeTopologyGroup
    (
        vtkhdfGroup_,
        "Lines",
        emptyConnectivity,
        emptyOffsets,
        0,
        emptyRows,
        HDF5::perPieceOffsetsRows(emptyRows)
    );
    writeTopologyGroup
    (
        vtkhdfGroup_,
        "Polygons",
        connectivity,
        offsets,
        static_cast<long long>(boundaryFaceIndices_.size()),
        connectivityRows,
        HDF5::perPieceOffsetsRows(polygonRows)
    );
    writeTopologyGroup
    (
        vtkhdfGroup_,
        "Strips",
        emptyConnectivity,
        emptyOffsets,
        0,
        emptyRows,
        HDF5::perPieceOffsetsRows(emptyRows)
    );

    globalNumBoundaryFaces_ = static_cast<Count>(polygonRows.globalRows);
    boundaryFaceRowOffset_ = static_cast<Index>(polygonRows.rowOffset);
    geometryWritten_ = true;

    closeHandles();

    if (debug_)
    {
        std::cout << std::format
        (
            "VTKHDF boundary geometry written: {}\n"
            "  - Boundary faces: {}\n"
            "  - Boundary nodes: {}\n"
            "  - Number of pieces: {}\n",
            fileName_,
            globalNumBoundaryFaces_,
            pointRows.globalRows,
            Comm::numProcessors()
        );
    }
}


void HDF5BoundaryData::appendTimeStep
(
    Scalar time,
    const FaceDataMap& scalarFaceFields
)
{
    if (!geometryWritten_)
    {
        FatalError
        (
            "VTKHDF appendTimeStep() called before writeGeometry(): "
          + fileName_
        );
    }

    if (finalized_)
    {
        FatalError
        (
            "VTKHDF appendTimeStep() called after finalize(): " + fileName_
        );
    }

    openHandles();

    const Count numBoundaryFaces = boundaryFaceIndices_.size();
    const FaceList& allFaces = mesh_.faces();

    // Static patch metadata, re-emitted every step
    checkFieldConsistency("patchIdx");
    appendFaceDataArray
    (
        "patchIdx",
        H5T_STD_I32LE,
        H5T_NATIVE_INT,
        patchIdx_.data()
    );
    checkFieldConsistency("patchZoneIdx");
    appendFaceDataArray
    (
        "patchZoneIdx",
        H5T_STD_I32LE,
        H5T_NATIVE_INT,
        patchZoneIdx_.data()
    );
    checkFieldConsistency("patchTypeIdx");
    appendFaceDataArray
    (
        "patchTypeIdx",
        H5T_STD_I32LE,
        H5T_NATIVE_INT,
        patchTypeIdx_.data()
    );
    checkFieldConsistency("isWall");
    appendFaceDataArray
    (
        "isWall",
        H5T_STD_U8LE,
        H5T_NATIVE_UCHAR,
        isWall_.data()
    );

    Count processedArrays = 4;

    // Per-step user face fields, reduced to boundary-face order
    std::vector<float> faceBuffer(numBoundaryFaces);

    for (const auto& [fieldName, faceField] : scalarFaceFields)
    {
        if (!faceField) continue;

        if (faceField->size() != allFaces.size())
        {
            FatalError
            (
                "VTKHDF boundary scalar field '" + fieldName + "' has "
              + std::to_string(faceField->size())
              + " values but the mesh has "
              + std::to_string(allFaces.size()) + " faces."
            );
        }

        checkFieldConsistency(fieldName);

        for (Index i = 0; i < numBoundaryFaces; ++i)
        {
            faceBuffer[i] =
                static_cast<float>((*faceField)[boundaryFaceIndices_[i]]);
        }

        appendFaceDataArray
        (
            fieldName,
            H5T_IEEE_F32LE,
            H5T_NATIVE_FLOAT,
            faceBuffer.data()
        );
        ++processedArrays;
    }

    if (stepCount_ > 0 && processedArrays != cellDataComponents_.size())
    {
        FatalError
        (
            "VTKHDF boundary field set changed between time steps: "
          + std::to_string(processedArrays) + " arrays written but "
          + std::to_string(cellDataComponents_.size())
          + " were registered on the first step."
        );
    }

    appendStepBookkeeping(time);
    ++stepCount_;

    // Refresh NSteps, then release the file until the next write
    HDF5::writeIntScalarAttribute
    (
        stepsGroup_,
        "NSteps",
        static_cast<long long>(stepCount_)
    );
    closeHandles();
}


void HDF5BoundaryData::finalize()
{
    if (finalized_)
    {
        return;
    }

    closeHandles();

    finalized_ = true;
}

// ****************************** Private Methods *****************************

void HDF5BoundaryData::openHandles()
{
    file_ = HDF5::openFile(fileName_);
    vtkhdfGroup_ = HDF5::openGroup(file_, "VTKHDF");
    cellDataGroup_ = HDF5::openGroup(vtkhdfGroup_, "CellData");
    stepsGroup_ = HDF5::openGroup(vtkhdfGroup_, "Steps");
    cellDataOffsetsGroup_ = HDF5::openGroup(stepsGroup_, "CellDataOffsets");
}


void HDF5BoundaryData::closeHandles()
{
    HDF5::closeGroup(cellDataOffsetsGroup_);
    HDF5::closeGroup(stepsGroup_);
    HDF5::closeGroup(cellDataGroup_);
    HDF5::closeGroup(vtkhdfGroup_);
    HDF5::closeFile(file_);
}


void HDF5BoundaryData::appendFaceDataArray
(
    const Name& name,
    hid_t fileType,
    hid_t memType,
    const void* values
)
{
    // Chunk rows derive from the GLOBAL count for the collective create
    const HDF5::SlabLayout dataRows
    {
        static_cast<hsize_t>(boundaryFaceIndices_.size()),
        static_cast<hsize_t>(globalNumBoundaryFaces_),
        static_cast<hsize_t>(boundaryFaceRowOffset_)
    };

    HDF5::appendRows
    (
        cellDataGroup_,
        name.c_str(),
        fileType,
        memType,
        values,
        dataRows,
        0,
        std::min<hsize_t>
        (
            std::max<hsize_t>(globalNumBoundaryFaces_, 1),
            HDF5::maxChunkRows
        )
    );

    const long long offsetValue =
        static_cast<long long>(cellDataRowOffsets_[name]);

    HDF5::appendRows
    (
        cellDataOffsetsGroup_,
        name.c_str(),
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &offsetValue,
        HDF5::replicatedRows(1),
        0,
        HDF5::stepChunkRows
    );

    cellDataRowOffsets_[name] += globalNumBoundaryFaces_;
}


void HDF5BoundaryData::checkFieldConsistency(const Name& name)
{
    if (stepCount_ == 0)
    {
        cellDataComponents_[name] = 1;
        cellDataRowOffsets_[name] = 0;
        return;
    }

    if (cellDataComponents_.find(name) == cellDataComponents_.end())
    {
        FatalError
        (
            "VTKHDF boundary field '" + name
          + "' appeared after the first time step; the field set must be "
            "identical on every step."
        );
    }
}


void HDF5BoundaryData::appendStepBookkeeping(Scalar time)
{
    const double timeValue = static_cast<double>(time);
    const long long zero = 0;
    const long long numParts = static_cast<long long>(Comm::numProcessors());
    const long long zero4[4] = {0, 0, 0, 0};

    // Identical on every rank: the master writes, the others just join
    const HDF5::SlabLayout oneRow = HDF5::replicatedRows(1);

    HDF5::appendRows
    (
        stepsGroup_,
        "Values",
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        &timeValue,
        oneRow,
        0,
        HDF5::stepChunkRows
    );

    // Static geometry: every step re-reads all parts from offset zero
    HDF5::appendRows
    (
        stepsGroup_,
        "PartOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &zero,
        oneRow,
        0,
        HDF5::stepChunkRows
    );
    HDF5::appendRows
    (
        stepsGroup_,
        "NumberOfParts",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numParts,
        oneRow,
        0,
        HDF5::stepChunkRows
    );
    HDF5::appendRows
    (
        stepsGroup_,
        "PointOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &zero,
        oneRow,
        0,
        HDF5::stepChunkRows
    );

    // [NSteps, 4]: Vertices, Lines, Polygons, Strips in that order
    HDF5::appendRows
    (
        stepsGroup_,
        "CellOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        zero4,
        oneRow,
        4,
        HDF5::stepChunkRows
    );
    HDF5::appendRows
    (
        stepsGroup_,
        "ConnectivityIdOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        zero4,
        oneRow,
        4,
        HDF5::stepChunkRows
    );
}

} // namespace VTK
