/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HDF5BoundaryData.h
 * @brief VTKHDF (.vtkhdf) boundary data for wall and patch face data
 *
 * @details Writes a static-geometry temporal VTKHDF PolyData holding all
 * boundary patch faces as polygons over a compressed local point set.
 *
 * Patch-metadata arrays are re-emitted every step because a temporal reader
 * expects an offsets entry for every cell-data array at every step. The file
 * is fully closed between writes.
 *
 * Parallel runs share ONE file: every rank holds a writer for the same
 * path and each rank's boundary faces (processor cuts excluded) form one
 * piece of a multi-piece layout, written collectively via MPI-IO. Every
 * method is collective — all ranks must call it, with an identical field
 * set, even when a rank's boundary is empty.
 *
 * @class HDF5BoundaryData
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <map>
#include <vector>

// External library headers
#include <hdf5.h>

// Project headers
#include "Scalar.h"
#include "Integer.h"
#include "StringTypes.h"
#include "Mesh.h"
#include "FaceData.h"

// *************************** Forward Declarations ***************************

class BoundaryConditions;

// ******************************* namespace VTK ******************************

namespace VTK
{

/// Face-centered scalar fields keyed by output name
using FaceDataMap = std::map<Name, const FaceData<Scalar>*>;

// ************************** class HDF5BoundaryData **************************

class HDF5BoundaryData
{
public:

// ************************* Special Member Functions *************************

    /// Create (truncate) the .vtkhdf file and its VTKHDF skeleton
    HDF5BoundaryData
    (
        FilePath fileName,
        const Mesh& mesh,
        const BoundaryConditions* bc = nullptr,
        bool debug = false
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    HDF5BoundaryData(const HDF5BoundaryData&) = delete;
    HDF5BoundaryData& operator=(const HDF5BoundaryData&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    HDF5BoundaryData(HDF5BoundaryData&&) = delete;
    HDF5BoundaryData& operator=(HDF5BoundaryData&&) = delete;

    /// Destructor - finalizes if not already done
    ~HDF5BoundaryData() noexcept;

// ****************************** Public Methods ******************************

    /// Write the static boundary polygon geometry once
    void writeGeometry();

    /// Append one time step of face-centered fields
    void appendTimeStep
    (
        Scalar time,
        const FaceDataMap& scalarFaceFields
    );

    /// Close all HDF5 handles
    void finalize();

    /// Number of time steps appended so far
    [[nodiscard]] Count stepCount() const noexcept
    {
        return stepCount_;
    }

// ****************************** Private Methods *****************************

private:

    /// Reopen the file and the VTKHDF groups for the next write
    void openHandles();

    /// Close all HDF5 handles
    void closeHandles();

    /// Append one cell-data array and its running offset row
    void appendFaceDataArray
    (
        const Name& name,
        hid_t fileType,
        hid_t memType,
        const void* values
    );

    /// Register (first step) or check (later steps) an array name
    void checkFieldConsistency(const Name& name);

    /// Append the per-step time value and all-zero geometry offset rows
    void appendStepBookkeeping(Scalar time);

// ****************************** Private Members *****************************

    /// Non-owning mesh reference
    const Mesh& mesh_;

    /// Optional boundary conditions view for wall identification
    const BoundaryConditions* bc_ = nullptr;

    /// Output file path
    FilePath fileName_;

    /// Verbose diagnostics flag
    bool debug_;

    /// HDF5 file handle (open only inside a write call)
    hid_t file_;

    /// /VTKHDF group
    hid_t vtkhdfGroup_;

    /// /VTKHDF/CellData group
    hid_t cellDataGroup_;

    /// /VTKHDF/Steps group
    hid_t stepsGroup_;

    /// /VTKHDF/Steps/CellDataOffsets group
    hid_t cellDataOffsetsGroup_;

    /// True once writeGeometry() has run
    bool geometryWritten_;

    /// True once finalize() has run
    bool finalized_;

    /// Number of time steps appended
    Count stepCount_;

    /// Global boundary-face count across every rank (set by writeGeometry)
    Count globalNumBoundaryFaces_;

    /// Global row where this rank's face slab starts (set by writeGeometry)
    Index boundaryFaceRowOffset_;

    /// Mesh face index of each exported boundary face, in polygon order
    IndexList boundaryFaceIndices_;

    /// Static per-face patch metadata, cached for per-step re-emission
    std::vector<int> patchIdx_;
    std::vector<int> patchZoneIdx_;
    std::vector<int> patchTypeIdx_;
    std::vector<unsigned char> isWall_;

    /// Per-array cumulative cell-data row offset
    std::map<Name, Count> cellDataRowOffsets_;

    /// Registered array names, for cross-step consistency checks
    std::map<Name, Count> cellDataComponents_;
};

} // namespace VTK
