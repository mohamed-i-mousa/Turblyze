#include "VtkWriter.h"
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cmath>

#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "Scalar.h"
#include "CellData.h"

namespace VtkWriter {

// Helper function to order nodes for a hexahedron
std::vector<size_t> orderHexNodes(const std::vector<size_t>& nodeIndices,
                                  const std::vector<Vector>& allNodes) 
    {
    // VTK_HEXAHEDRON expects the following node ordering (see VTK file format):
    //   0 (-x,-y,-z)
    //   1 (+x,-y,-z)
    //   2 (+x,+y,-z)
    //   3 (-x,+y,-z)
    //   4 (-x,-y,+z)
    //   5 (+x,-y,+z)
    //   6 (+x,+y,+z)
    //   7 (-x,+y,+z)

    constexpr double TOLERANCE = 1e-8;

    // --- Separate bottom (min-z) and top (max-z) layers --- //
    double minZ = std::numeric_limits<double>::max();
    double maxZ = std::numeric_limits<double>::lowest();
    for (size_t idx : nodeIndices) {
        double z = allNodes[idx].z;
        minZ = std::min(minZ, z);
        maxZ = std::max(maxZ, z);
    }

    std::vector<size_t> bottomLayer;
    std::vector<size_t> topLayer;
    for (size_t idx : nodeIndices) {
        double z = allNodes[idx].z;
        if (std::abs(z - minZ) < TOLERANCE) {
            bottomLayer.push_back(idx);
        } else if (std::abs(z - maxZ) < TOLERANCE) {
            topLayer.push_back(idx);
        }
    }

    // Sanity check: fall back to original ordering if something goes wrong
    if (bottomLayer.size() != 4 || topLayer.size() != 4) {
        return nodeIndices;
    }

    // --- Order bottom layer counter-clockwise around its centroid --- //
    Vector cBottom(0, 0, 0);
    for (size_t idx : bottomLayer) cBottom += allNodes[idx];
    cBottom = cBottom * (1.0 / bottomLayer.size());

    auto angleCmp = [&](size_t a, size_t b) {
        double angA = std::atan2(allNodes[a].y - cBottom.y, allNodes[a].x - cBottom.x);
        double angB = std::atan2(allNodes[b].y - cBottom.y, allNodes[b].x - cBottom.x);
        return angA < angB;
    };
    std::sort(bottomLayer.begin(), bottomLayer.end(), angleCmp);

    // Rotate so that the first bottom node is the one with smallest (x,y)
    auto firstIt = std::min_element(bottomLayer.begin(), bottomLayer.end(), [&](size_t a, size_t b) {
        const Vector& pa = allNodes[a];
        const Vector& pb = allNodes[b];
        if (pa.x != pb.x) return pa.x < pb.x;
        return pa.y < pb.y;
    });
    std::rotate(bottomLayer.begin(), firstIt, bottomLayer.end());

    // --- Order the top layer so that it matches the bottom ordering --- //
    std::vector<size_t> topOrdered;
    topOrdered.reserve(4);
    for (size_t bIdx : bottomLayer) {
        const Vector& b = allNodes[bIdx];
        size_t bestIdx = topLayer[0];
        double bestDist2 = std::numeric_limits<double>::max();
        for (size_t tIdx : topLayer) {
            double dx = allNodes[tIdx].x - b.x;
            double dy = allNodes[tIdx].y - b.y;
            double dist2 = dx * dx + dy * dy;
            if (dist2 < bestDist2) {
                bestDist2 = dist2;
                bestIdx = tIdx;
            }
        }
        topOrdered.push_back(bestIdx);
    }

    // --- Combine layers in VTK order --- //
    std::vector<size_t> result;
    result.reserve(8);
    result.insert(result.end(), bottomLayer.begin(), bottomLayer.end());
    result.insert(result.end(), topOrdered.begin(), topOrdered.end());
    return result;
}

// Helper function to order nodes for a tetrahedron
std::vector<size_t> orderTetNodes(const std::vector<size_t>& nodeIndices,
                                  const std::vector<Vector>& allNodes) {
    // Find centroid
    Vector centroid(0, 0, 0);
    for (size_t nodeIdx : nodeIndices) {
        centroid = centroid + allNodes[nodeIdx];
    }
    centroid = centroid * (1.0 / nodeIndices.size());

    // Find the base node (farthest from centroid)
    size_t baseNode = nodeIndices[0];
    double maxDist = (allNodes[baseNode] - centroid).magnitudeSquared();
    for (size_t i = 1; i < nodeIndices.size(); i++) {
        double dist = (allNodes[nodeIndices[i]] - centroid).magnitudeSquared();
        if (dist > maxDist) {
            maxDist = dist;
            baseNode = nodeIndices[i];
        }
    }

    // Order the other three nodes to form the base
    std::vector<size_t> otherNodes;
    for (size_t nodeIdx : nodeIndices) {
        if (nodeIdx != baseNode) {
            otherNodes.push_back(nodeIdx);
        }
    }

    // Ensure counter-clockwise order for the base
    Vector v1 = allNodes[otherNodes[1]] - allNodes[otherNodes[0]];
    Vector v2 = allNodes[otherNodes[2]] - allNodes[otherNodes[0]];
    Vector normal = cross(v1, v2);
    
    // Check if the order needs to be reversed
    Vector toCentroid = centroid - allNodes[otherNodes[0]];
    if (dot(normal, toCentroid) < 0) {
        std::swap(otherNodes[1], otherNodes[2]);
    }

    // VTK order: base nodes first, then apex
    return {otherNodes[0], otherNodes[1], otherNodes[2], baseNode};
}

// Helper function to order nodes for a pyramid (square base + apex)
std::vector<size_t> orderPyramidNodes(const std::vector<size_t>& nodeIndices,
                                     const std::vector<Vector>& allNodes) {
    // Identify apex: the node farthest from the average of all nodes
    Vector centroid(0, 0, 0);
    for (size_t idx : nodeIndices) centroid += allNodes[idx];
    centroid = centroid * (1.0 / nodeIndices.size());

    size_t apexIdx = nodeIndices[0];
    double maxDist = (allNodes[apexIdx] - centroid).magnitudeSquared();
    for (size_t idx : nodeIndices) {
        double dist = (allNodes[idx] - centroid).magnitudeSquared();
        if (dist > maxDist) { apexIdx = idx; maxDist = dist; }
    }

    // Base nodes are the remaining 4
    std::vector<size_t> baseNodes;
    for (size_t idx : nodeIndices) if (idx != apexIdx) baseNodes.push_back(idx);

    // Order base nodes counter-clockwise around their centroid
    Vector baseCentroid(0,0,0);
    for (size_t idx : baseNodes) baseCentroid += allNodes[idx];
    baseCentroid = baseCentroid * (1.0 / baseNodes.size());

    auto angleCmp = [&](size_t a, size_t b) {
        double angA = std::atan2(allNodes[a].y - baseCentroid.y,
                                 allNodes[a].x - baseCentroid.x);
        double angB = std::atan2(allNodes[b].y - baseCentroid.y,
                                 allNodes[b].x - baseCentroid.x);
        return angA < angB;
    };
    std::sort(baseNodes.begin(), baseNodes.end(), angleCmp);

    // VTK order: base square (counter-clockwise) then apex
    std::vector<size_t> result = baseNodes;
    result.push_back(apexIdx);
    return result;
}

// Helper function to order nodes for a wedge (triangular prism)
std::vector<size_t> orderWedgeNodes(const std::vector<size_t>& nodeIndices,
                                    const std::vector<Vector>& allNodes) {
    // Split nodes into two layers based on Z (or distance along principal axis)
    constexpr double TOL = 1e-8;
    double minZ = std::numeric_limits<double>::max();
    double maxZ = std::numeric_limits<double>::lowest();
    for (size_t idx : nodeIndices) {
        double z = allNodes[idx].z;
        minZ = std::min(minZ, z);
        maxZ = std::max(maxZ, z);
    }
    std::vector<size_t> lower, upper;
    for (size_t idx : nodeIndices) {
        double z = allNodes[idx].z;
        if (std::abs(z - minZ) < TOL) lower.push_back(idx); else upper.push_back(idx);
    }
    // Fallback if not split correctly
    if (lower.size() != 3 || upper.size() != 3) return nodeIndices;

    // Order the lower triangle counter-clockwise
    Vector cLower(0,0,0);
    for (size_t idx : lower) cLower += allNodes[idx];
    cLower = cLower * (1.0 / lower.size());
    auto angCmp = [&](size_t a, size_t b) {
        double angA = std::atan2(allNodes[a].y - cLower.y, allNodes[a].x - cLower.x);
        double angB = std::atan2(allNodes[b].y - cLower.y, allNodes[b].x - cLower.x);
        return angA < angB;
    };
    std::sort(lower.begin(), lower.end(), angCmp);

    // Match upper nodes to lower ones (closest in x-y)
    std::vector<size_t> upperOrdered;
    upperOrdered.reserve(3);
    for (size_t idxL : lower) {
        size_t best = upper[0]; double bestDist2 = std::numeric_limits<double>::max();
        for (size_t idxU : upper) {
            double dx = allNodes[idxU].x - allNodes[idxL].x;
            double dy = allNodes[idxU].y - allNodes[idxL].y;
            double d2 = dx*dx + dy*dy;
            if (d2 < bestDist2) { bestDist2 = d2; best = idxU; }
        }
        upperOrdered.push_back(best);
    }

    // VTK order: lower triangle (0-2), upper triangle (3-5)
    std::vector<size_t> result;
    result.reserve(6);
    result.insert(result.end(), lower.begin(), lower.end());
    result.insert(result.end(), upperOrdered.begin(), upperOrdered.end());
    return result;
}

void writeVtkFile(const std::string& filename,
                 const std::vector<Vector>& allNodes,
                 const std::vector<Face>& allFaces,
                 const std::vector<Cell>& allCells,
                 const std::map<std::string, const ScalarField*>& scalarCellFields) {
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open()) {
        throw std::runtime_error("Error: Could not open VTK file for writing: " + filename);
    }

    // --- VTK Header ---
    vtkFile << std::fixed << std::setprecision(8);
    vtkFile << "# vtk DataFile Version 4.2\n";
    vtkFile << "Unstructured Grid Data\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";

    // --- Points Section ---
    vtkFile << "POINTS " << allNodes.size() << " ";
    if constexpr (std::is_same_v<Scalar, double>) {
        vtkFile << "double\n";
    } else {
        vtkFile << "float\n";
    }
    for (const auto& node : allNodes) {
        vtkFile << node.x << " " << node.y << " " << node.z << "\n";
    }
    vtkFile << "\n";

    // --- Cells Section ---
    std::vector<const Cell*> validCells;
    for (const auto& cell : allCells) {
        if (!cell.faceIndices.empty()) {
            validCells.push_back(&cell);
        }
    }

    const size_t numCells = validCells.size();

    // Calculate total data size for CELLS section and validate indices
    size_t totalCellDataSize = 0;
    std::vector<int> cellTypes(numCells, 42);  // Default to polyhedron (42)
    std::vector<std::vector<size_t>> orderedNodes(numCells);

    for (size_t i = 0; i < numCells; i++) {
        const Cell& cell = *validCells[i];
        std::unordered_set<size_t> uniqueNodes;

        // Collect unique nodes and validate indices
        for (size_t faceIdx : cell.faceIndices) {
            if (faceIdx >= allFaces.size()) {
                throw std::out_of_range("Face index " + std::to_string(faceIdx) +
                                       " in cell " + std::to_string(cell.id) +
                                       " is out of bounds");
            }
            const Face& face = allFaces[faceIdx];
            for (size_t nodeIdx : face.nodeIndices) {
                if (nodeIdx >= allNodes.size()) {
                    throw std::out_of_range("Node index " + std::to_string(nodeIdx) +
                                           " in cell " + std::to_string(cell.id) +
                                           " exceeds points list size");
                }
                uniqueNodes.insert(nodeIdx);
            }
        }

        // Determine cell type based on topology
        const size_t numNodes = uniqueNodes.size();
        const size_t numFaces = cell.faceIndices.size();

        // Create vector of unique node indices
        std::vector<size_t> nodeIndices(uniqueNodes.begin(), uniqueNodes.end());

        if (numNodes == 4 && numFaces == 4) {
            // Check if all faces are triangles
            bool allTriangles = true;
            for (size_t faceIdx : cell.faceIndices) {
                if (allFaces[faceIdx].nodeIndices.size() != 3) {
                    allTriangles = false;
                    break;
                }
            }
            if (allTriangles) {
                cellTypes[i] = 10;  // VTK_TETRA
                orderedNodes[i] = orderTetNodes(nodeIndices, allNodes);
                totalCellDataSize += 1 + 4;  // 1 (size) + 4 (nodes)
                continue;
            }
        } else if (numNodes == 8 && numFaces == 6) {
            // Check if all faces are quadrilaterals
            bool allQuads = true;
            for (size_t faceIdx : cell.faceIndices) {
                if (allFaces[faceIdx].nodeIndices.size() != 4) {
                    allQuads = false;
                    break;
                }
            }
            if (allQuads) {
                cellTypes[i] = 12;  // VTK_HEXAHEDRON
                orderedNodes[i] = orderHexNodes(nodeIndices, allNodes);
                totalCellDataSize += 1 + 8;  // 1 (size) + 8 (nodes)
                continue;
            }
        } else if (numNodes == 5 && numFaces == 5) {
            // Likely a pyramid: 1 quad base + 4 tri faces
            bool hasQuadFace = false;
            for (size_t faceIdx : cell.faceIndices) {
                if (allFaces[faceIdx].nodeIndices.size() == 4) {
                    hasQuadFace = true;
                    break;
                }
            }
            if (hasQuadFace) {
                cellTypes[i] = 14;  // VTK_PYRAMID
                orderedNodes[i] = orderPyramidNodes(nodeIndices, allNodes);
                totalCellDataSize += 1 + 5;  // 1 (size) + 5 (nodes)
                continue;
            }
        } else if (numNodes == 6 && numFaces == 5) {
            // Likely a wedge (triangular prism): 2 tri & 3 quad faces
            cellTypes[i] = 13;  // VTK_WEDGE
            orderedNodes[i] = orderWedgeNodes(nodeIndices, allNodes);
            totalCellDataSize += 1 + 6;  // 1 (size) + 6 (nodes)
            continue;
        }

        // Polyhedral cell
        cellTypes[i] = 42;
        size_t faceDataSize = 1;  // For numFaces
        for (size_t faceIdx : cell.faceIndices) {
            const Face& face = allFaces[faceIdx];
            faceDataSize += 1 + face.nodeIndices.size();  // Face size + nodes
        }
        totalCellDataSize += 1 + faceDataSize;  // Cell size + face data
    }

    // Write CELLS section
    vtkFile << "CELLS " << numCells << " " << totalCellDataSize << "\n";

    for (size_t i = 0; i < numCells; i++) {
        const Cell& cell = *validCells[i];
        const int cellType = cellTypes[i];

        if (cellType == 10) {  // Tetrahedron
            vtkFile << "4";
            for (auto nodeIdx : orderedNodes[i]) {
                vtkFile << " " << nodeIdx;
            }
            vtkFile << "\n";
        } else if (cellType == 12) {  // Hexahedron
            vtkFile << "8";
            for (auto nodeIdx : orderedNodes[i]) {
                vtkFile << " " << nodeIdx;
            }
            vtkFile << "\n";
        } else if (cellType == 14) {  // Pyramid
            vtkFile << "5";
            for (auto nodeIdx : orderedNodes[i]) {
                vtkFile << " " << nodeIdx;
            }
            vtkFile << "\n";
        } else if (cellType == 13) {  // Wedge
            vtkFile << "6";
            for (auto nodeIdx : orderedNodes[i]) {
                vtkFile << " " << nodeIdx;
            }
            vtkFile << "\n";
        } else {  // Polyhedron
            vtkFile << cell.faceIndices.size();  // Number of faces

            for (size_t faceIdx : cell.faceIndices) {
                const Face& face = allFaces[faceIdx];
                vtkFile << " " << face.nodeIndices.size();  // Points in face
                for (auto nodeIdx : face.nodeIndices) {
                    vtkFile << " " << nodeIdx;
                }
            }
            vtkFile << "\n";
        }
    }

    // ----- Cell Types Section ----- //
    vtkFile << "CELL_TYPES " << numCells << "\n";
    for (size_t i = 0; i < numCells; ++i) {
        vtkFile << cellTypes[i] << "\n";
    }

    // ----- Cell Data Section ----- //
    if (numCells > 0 && !scalarCellFields.empty()) {
        vtkFile << "CELL_DATA " << numCells << "\n";

        // Write Scalar Fields
        for (const auto& pair : scalarCellFields) {
            const std::string& fieldName = pair.first;
            const ScalarField* field = pair.second;

            if (field && field->size() == allCells.size()) {
                vtkFile << "SCALARS " << fieldName << " ";
                if constexpr (std::is_same_v<Scalar, double>) {
                    vtkFile << "double 1\n";
                } else {
                    vtkFile << "float 1\n";
                }
                vtkFile << "LOOKUP_TABLE default\n";
                for (const auto* cellPtr : validCells) {
                    vtkFile << (*field)[cellPtr->id] << "\n";
                }
                vtkFile << "\n";
            } else {
                std::cerr << "Warning (VTK): Scalar field '" << fieldName 
                          << "' is null or size mismatch. Skipping." << std::endl;
            }
        }
    }

    vtkFile.close();
    if (vtkFile.fail()) {
        throw std::runtime_error("Error closing file: " + filename);
    }
}

} // namespace VtkWriter