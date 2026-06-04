/******************************************************************************
 * @file MeshCreator.cpp
 * @brief Mesh read and geometry-preparation phase
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "MeshCreator.h"

// Standard library headers
#include <iostream>
#include <vector>

// Project headers
#include "Logger.h"
#include "MeshChecker.h"
#include "MeshReader.h"

// *************************** namespace MeshCreator **************************

namespace MeshCreator
{

Mesh create(const CaseConfiguration& config)
{
    std::cout << '\n';
    Logger::sectionHeader("Reading and Preparing Mesh");

    // Read mesh data from file
    MeshReader meshReader(config.meshFile);

    // Construct Mesh
    Mesh mesh
    (
        meshReader.moveNodes(),
        meshReader.moveFaces(),
        meshReader.moveCells(),
        meshReader.moveBoundaryPatches()
    );

    std::cout
        << "Mesh Loaded: " << mesh.numNodes() << " nodes, "
        << mesh.numFaces() << " faces, " << mesh.numCells()
        << " cells." << '\n';

    // Calculate geometric properties for faces
    std::vector<FaceIntegrals> faceIntegrals(mesh.numFaces());
    auto faces = mesh.faces();
    const auto nodes = mesh.nodes();

    #pragma omp parallel for schedule(static)
    for (Index faceIdx = 0; faceIdx < faces.size(); ++faceIdx)
    {
        faceIntegrals[faceIdx] =
            faces[faceIdx].geometricProperties(nodes);
    }
    if (config.debug)
    {
        std::cout
            << "Geometric properties calculated for faces." << '\n';
    }

    // Calculate geometric properties for cells
    auto cells = mesh.cells();

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < cells.size(); ++cellIdx)
    {
        cells[cellIdx].geometricProperties(faces, faceIntegrals);
    }
    if (config.debug)
    {
        std::cout
            << "Geometric properties calculated for cells." << '\n';
    }

    #pragma omp parallel for schedule(static)
    for (Index faceIdx = 0; faceIdx < faces.size(); ++faceIdx)
    {
        faces[faceIdx].distances(cells);
    }
    if (config.debug)
    {
        std::cout
            << "Distance properties calculated for faces." << '\n';
    }

    // Check mesh quality if enabled
    if (config.checkQuality)
    {
        MeshChecker::check(mesh);
    }

    return mesh;
}

} // namespace MeshCreator
