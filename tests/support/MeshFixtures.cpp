/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshFixtures.cpp
 * @brief Implementation of the programmatic hex-box mesh builder
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "MeshFixtures.h"

// Standard library headers
#include <cstdint>
#include <utility>
#include <vector>

// Project headers
#include "Cell.h"
#include "Face.h"
#include "BoundaryPatch.h"
#include "MeshContainers.h"
#include "ProcessorPatch.h"
#include "MeshCreator.h"
#include "Comm.h"
#include "HaloExchange.h"

// ***************************** Internal Helpers *****************************

namespace
{

/// Per-cell accumulators
struct CellTopology
{
    std::vector<IndexList> faces;
    std::vector<IndexList> neighbors;
    std::vector<Cell::FaceSignList> signs;
};

/// Linear index of node (i, j, k) in an (nx+1)(ny+1)(nz+1) grid
[[nodiscard]] Index nodeIndex
(
    Index i,
    Index j,
    Index k,
    Count nx,
    Count ny
) noexcept
{
    return i + (nx + 1) * (j + (ny + 1) * k);
}

/// Linear index of cell (i, j, k) in an nx*ny*nz grid
[[nodiscard]] Index cellIndex
(
    Index i,
    Index j,
    Index k,
    Count nx,
    Count ny
) noexcept
{
    return i + nx * (j + ny * k);
}

/// Append an internal face and register it on both cells
void addInternalFace
(
    FaceList& faces,
    CellTopology& topology,
    IndexList nodeIds,
    Index owner,
    Index neighbor
)
{
    const Index faceIdx = faces.size();

    faces.emplace_back(faceIdx, std::move(nodeIds), owner, neighbor);

    topology.faces[owner].push_back(faceIdx);
    topology.signs[owner].push_back(static_cast<std::int8_t>(1));
    topology.neighbors[owner].push_back(neighbor);

    topology.faces[neighbor].push_back(faceIdx);
    topology.signs[neighbor].push_back(static_cast<std::int8_t>(-1));
    topology.neighbors[neighbor].push_back(owner);
}

/// Append a boundary face and register it on its owner with sign +1
void addBoundaryFace
(
    FaceList& faces,
    CellTopology& topology,
    IndexList nodeIds,
    Index owner
)
{
    const Index faceIdx = faces.size();

    faces.emplace_back(faceIdx, std::move(nodeIds), owner);

    topology.faces[owner].push_back(faceIdx);
    topology.signs[owner].push_back(static_cast<std::int8_t>(1));
}

/// Append a boundary patch
void addBoundaryPatch
(
    PatchList& patches,
    Index zoneIdx,
    Index firstFaceIdx,
    const FaceList& faces,
    const Name& patchName
)
{
    BoundaryPatch patch(zoneIdx, firstFaceIdx, faces.size() - 1);
    patch.setName(patchName);
    patch.setType(PatchType::physical);
    patches.push_back(std::move(patch));
}

} // namespace

// **************************** Fixture Factories *****************************

Mesh makeHexBoxMesh
(
    Count nx,
    Count ny,
    Count nz,
    Scalar spacing
)
{
    // Nodes: grid points pushed in ascending linear-index order
    NodeList nodes;
    nodes.reserve((nx + 1) * (ny + 1) * (nz + 1));

    for (Index k = 0; k <= nz; ++k)
    {
        for (Index j = 0; j <= ny; ++j)
        {
            for (Index i = 0; i <= nx; ++i)
            {
                nodes.emplace_back
                (
                    S(i) * spacing,
                    S(j) * spacing,
                    S(k) * spacing
                );
            }
        }
    }

    const Count numCells = nx * ny * nz;

    CellTopology topology;
    topology.faces.resize(numCells);
    topology.neighbors.resize(numCells);
    topology.signs.resize(numCells);

    FaceList faces;

    // Internal faces first (so boundary patches stay contiguous later)

    // x-normal internal faces at planes i = 1 .. nx-1 (+x winding, out of low)
    for (Index k = 0; k < nz; ++k)
    {
        for (Index j = 0; j < ny; ++j)
        {
            for (Index i = 1; i < nx; ++i)
            {
                addInternalFace
                (
                    faces,
                    topology,
                    IndexList
                    {
                        nodeIndex(i, j,     k,     nx, ny),
                        nodeIndex(i, j,     k + 1, nx, ny),
                        nodeIndex(i, j + 1, k + 1, nx, ny),
                        nodeIndex(i, j + 1, k,     nx, ny)
                    },
                    cellIndex(i - 1, j, k, nx, ny),
                    cellIndex(i,     j, k, nx, ny)
                );
            }
        }
    }

    // y-normal internal faces at planes j = 1 .. ny-1 (+y winding, out of low)
    for (Index k = 0; k < nz; ++k)
    {
        for (Index j = 1; j < ny; ++j)
        {
            for (Index i = 0; i < nx; ++i)
            {
                addInternalFace
                (
                    faces,
                    topology,
                    IndexList
                    {
                        nodeIndex(i,     j, k,     nx, ny),
                        nodeIndex(i + 1, j, k,     nx, ny),
                        nodeIndex(i + 1, j, k + 1, nx, ny),
                        nodeIndex(i,     j, k + 1, nx, ny)
                    },
                    cellIndex(i, j - 1, k, nx, ny),
                    cellIndex(i, j,     k, nx, ny)
                );
            }
        }
    }

    // z-normal internal faces at planes k = 1 .. nz-1 (+z winding, out of low)
    for (Index k = 1; k < nz; ++k)
    {
        for (Index j = 0; j < ny; ++j)
        {
            for (Index i = 0; i < nx; ++i)
            {
                addInternalFace
                (
                    faces,
                    topology,
                    IndexList
                    {
                        nodeIndex(i,     j,     k, nx, ny),
                        nodeIndex(i,     j + 1, k, nx, ny),
                        nodeIndex(i + 1, j + 1, k, nx, ny),
                        nodeIndex(i + 1, j,     k, nx, ny)
                    },
                    cellIndex(i, j, k - 1, nx, ny),
                    cellIndex(i, j, k,     nx, ny)
                );
            }
        }
    }

    // Boundary faces, one contiguous range per patch
    PatchList patches;

    // xMin (i = 0): outward normal -x, so the +x winding is reversed
    Index first = faces.size();
    for (Index k = 0; k < nz; ++k)
    {
        for (Index j = 0; j < ny; ++j)
        {
            addBoundaryFace
            (
                faces,
                topology,
                IndexList
                {
                    nodeIndex(0, j + 1, k,     nx, ny),
                    nodeIndex(0, j + 1, k + 1, nx, ny),
                    nodeIndex(0, j,     k + 1, nx, ny),
                    nodeIndex(0, j,     k,     nx, ny)
                },
                cellIndex(0, j, k, nx, ny)
            );
        }
    }
    addBoundaryPatch(patches, 0, first, faces, BoxPatch::xMin);

    // xMax (i = nx): outward normal +x
    first = faces.size();
    for (Index k = 0; k < nz; ++k)
    {
        for (Index j = 0; j < ny; ++j)
        {
            addBoundaryFace
            (
                faces,
                topology,
                IndexList
                {
                    nodeIndex(nx, j,     k,     nx, ny),
                    nodeIndex(nx, j,     k + 1, nx, ny),
                    nodeIndex(nx, j + 1, k + 1, nx, ny),
                    nodeIndex(nx, j + 1, k,     nx, ny)
                },
                cellIndex(nx - 1, j, k, nx, ny)
            );
        }
    }
    addBoundaryPatch(patches, 1, first, faces, BoxPatch::xMax);

    // yMin (j = 0): outward normal -y, reversed +y winding
    first = faces.size();
    for (Index k = 0; k < nz; ++k)
    {
        for (Index i = 0; i < nx; ++i)
        {
            addBoundaryFace
            (
                faces,
                topology,
                IndexList
                {
                    nodeIndex(i,     0, k + 1, nx, ny),
                    nodeIndex(i + 1, 0, k + 1, nx, ny),
                    nodeIndex(i + 1, 0, k,     nx, ny),
                    nodeIndex(i,     0, k,     nx, ny)
                },
                cellIndex(i, 0, k, nx, ny)
            );
        }
    }
    addBoundaryPatch(patches, 2, first, faces, BoxPatch::yMin);

    // yMax (j = ny): outward normal +y
    first = faces.size();
    for (Index k = 0; k < nz; ++k)
    {
        for (Index i = 0; i < nx; ++i)
        {
            addBoundaryFace
            (
                faces,
                topology,
                IndexList
                {
                    nodeIndex(i,     ny, k,     nx, ny),
                    nodeIndex(i + 1, ny, k,     nx, ny),
                    nodeIndex(i + 1, ny, k + 1, nx, ny),
                    nodeIndex(i,     ny, k + 1, nx, ny)
                },
                cellIndex(i, ny - 1, k, nx, ny)
            );
        }
    }
    addBoundaryPatch(patches, 3, first, faces, BoxPatch::yMax);

    // zMin (k = 0): outward normal -z, reversed +z winding
    first = faces.size();
    for (Index j = 0; j < ny; ++j)
    {
        for (Index i = 0; i < nx; ++i)
        {
            addBoundaryFace
            (
                faces,
                topology,
                IndexList
                {
                    nodeIndex(i + 1, j,     0, nx, ny),
                    nodeIndex(i + 1, j + 1, 0, nx, ny),
                    nodeIndex(i,     j + 1, 0, nx, ny),
                    nodeIndex(i,     j,     0, nx, ny)
                },
                cellIndex(i, j, 0, nx, ny)
            );
        }
    }
    addBoundaryPatch(patches, 4, first, faces, BoxPatch::zMin);

    // zMax (k = nz): outward normal +z
    first = faces.size();
    for (Index j = 0; j < ny; ++j)
    {
        for (Index i = 0; i < nx; ++i)
        {
            addBoundaryFace
            (
                faces,
                topology,
                IndexList
                {
                    nodeIndex(i,     j,     nz, nx, ny),
                    nodeIndex(i,     j + 1, nz, nx, ny),
                    nodeIndex(i + 1, j + 1, nz, nx, ny),
                    nodeIndex(i + 1, j,     nz, nx, ny)
                },
                cellIndex(i, j, nz - 1, nx, ny)
            );
        }
    }
    addBoundaryPatch(patches, 5, first, faces, BoxPatch::zMax);

    // Cells from the accumulated topology
    CellList cells;
    cells.reserve(numCells);

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        cells.emplace_back
        (
            cellIdx,
            std::move(topology.faces[cellIdx]),
            std::move(topology.neighbors[cellIdx]),
            std::move(topology.signs[cellIdx])
        );
    }

    MeshCreator::linkBoundaryFaces(faces, patches);

    // Geometry pass, mirroring MeshCreator::prepareGeometry
    std::vector<FaceIntegrals> faceIntegrals(faces.size());

    for (Index faceIdx = 0; faceIdx < faces.size(); ++faceIdx)
    {
        faceIntegrals[faceIdx] =
            faces[faceIdx].geometricProperties(nodes);
    }

    for (Index cellIdx = 0; cellIdx < cells.size(); ++cellIdx)
    {
        cells[cellIdx].geometricProperties(faceIntegrals);
    }

    return Mesh
    (
        std::move(nodes),
        std::move(faces),
        std::move(cells),
        std::move(patches)
    );
}


Mesh makeDecomposedChainMesh()
{
    const Index rank = Comm::myProcessorNum();
    const Count size = Comm::numProcessors();

    // Neighbours in ascending rank order, each becomes one ghost cell
    IndexList neighbors;
    if (rank > 0)
    {
        neighbors.push_back(rank - 1);
    }
    if (rank + 1 < size)
    {
        neighbors.push_back(rank + 1);
    }

    const Count numOwned = 1;
    const Count numGhost = neighbors.size();

    // One unit-square cut face per neighbour, fanned out from owned cell 0
    NodeList nodes;
    nodes.reserve(4 * numGhost);

    FaceList faces;
    faces.reserve(numGhost);

    IndexList ownedFaces;
    IndexList ownedNeighbors;
    Cell::FaceSignList ownedSigns;

    for (Index g = 0; g < numGhost; ++g)
    {
        const Scalar x = S(1 + 2 * g);
        const Index firstNode = nodes.size();

        nodes.emplace_back(x, S(0), S(0));
        nodes.emplace_back(x, S(1), S(0));
        nodes.emplace_back(x, S(1), S(1));
        nodes.emplace_back(x, S(0), S(1));

        const Index faceIdx = faces.size();

        faces.emplace_back
        (
            faceIdx,
            IndexList{firstNode, firstNode + 1, firstNode + 2, firstNode + 3},
            Index{0},               // owner: the owned cell
            numOwned + g            // neighbour: the ghost stub
        );

        ownedFaces.push_back(faceIdx);
        ownedNeighbors.push_back(numOwned + g);
        ownedSigns.push_back(static_cast<std::int8_t>(1));
    }

    // The owned cell carries the cut faces; the ghosts are geometry-only stubs
    CellList cells;
    cells.reserve(numOwned + numGhost);

    cells.emplace_back
    (
        Index{0},
        std::move(ownedFaces),
        std::move(ownedNeighbors),
        std::move(ownedSigns)
    );
    cells[0].setGeometry(Vector(S(0), S(0.5), S(0.5)), S(1));

    for (Index g = 0; g < numGhost; ++g)
    {
        Cell ghost;
        ghost.setIdx(numOwned + g);
        ghost.setGeometry(Vector(S(2 + 2 * g), S(0.5), S(0.5)), S(1));
        cells.push_back(std::move(ghost));
    }

    // One processor patch per neighbour
    PatchList patches;
    patches.reserve(numGhost);

    ProcessorPatchList processorPatches;
    processorPatches.reserve(numGhost);

    for (Index g = 0; g < numGhost; ++g)
    {
        BoundaryPatch patch(0, g, g);
        patch.setName
        (
            "processor" + std::to_string(rank)
          + "to" + std::to_string(neighbors[g])
        );
        patch.setType(PatchType::processor);
        patches.push_back(std::move(patch));

        processorPatches.emplace_back
        (
            neighbors[g],           // neighbour rank
            g,                      // first cut face
            g,                      // last cut face
            IndexList{0},           // send owned cell 0
            numOwned + g,           // ghost-tail slot this patch fills
            Count{1}                // one ghost cell
        );
    }

    MeshCreator::linkBoundaryFaces(faces, patches);

    // Face geometry only; the cells carry the centroids set above
    for (Index faceIdx = 0; faceIdx < faces.size(); ++faceIdx)
    {
        Face& face = faces[faceIdx];

        static_cast<void>(face.geometricProperties(nodes));
    }

    Mesh mesh
    (
        std::move(nodes),
        std::move(faces),
        std::move(cells),
        std::move(patches),
        numGhost
    );

    Halo::init
    (
        std::move(processorPatches),
        neighbors
    );

    return mesh;
}


Mesh makeDecomposedHexBoxMesh(Count nx, Count ny, Count nz, Scalar spacing)
{
    if (!Comm::parallelRun())
    {
        Halo::reset();
        return makeHexBoxMesh(nx, ny, nz, spacing);
    }

    // The production decomposition path, driven from the programmatic box
    Mesh completeMesh =
        Comm::master() ? makeHexBoxMesh(nx, ny, nz, spacing) : Mesh();

    return MeshCreator::decomposeAndDistribute(std::move(completeMesh), false);
}