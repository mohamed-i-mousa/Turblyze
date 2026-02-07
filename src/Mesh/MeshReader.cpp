/******************************************************************************
 * @file MeshReader.cpp
 * @brief Implementation of the MeshReader class for Fluent mesh files
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <charconv>

#include "MeshReader.hpp"


// ************************ Static member definitions *************************

const std::string MeshReader::MSH_COMMENT    = "(0";
const std::string MeshReader::MSH_DIMENSION  = "(2";
const std::string MeshReader::MSH_NODES      = "(10";
const std::string MeshReader::MSH_CELLS      = "(12";
const std::string MeshReader::MSH_FACES      = "(13";
const std::string MeshReader::MSH_BOUNDARIES = "(45";


// ******************************* Constructor ********************************

MeshReader::MeshReader(const std::string& filePath)
{
    parseFile(filePath);
}


// ************************* Static utility methods ***************************

size_t MeshReader::hexToDec(const std::string &hexStr)
{
    size_t decVal;
    auto result =
        std::from_chars
            (hexStr.data(), hexStr.data() + hexStr.size(), decVal, 16);

    if
    (
        result.ec != std::errc()
     || result.ptr != hexStr.data() + hexStr.size()
    )
    {
        throw   std::runtime_error
                (
                    "Error: Failed to convert hex string '" + hexStr
                  + "' to decimal."
                );
    }

    return decVal;
}

size_t MeshReader::strToDec(const std::string& decStr)
{
    if (decStr.empty())
    {
        throw   std::runtime_error
                (
                    "Error: Attempted to convert empty decimal string to "
                    "size_t."
                );
    }

    std::stringstream ss;
    ss  << decStr;
    size_t decVal;
    ss  >> decVal;

    if (ss.fail() || !ss.eof())
    {
        throw   std::runtime_error
                (
                    "Error: Failed to convert decimal string '" + decStr
                  + "' to size_t. Invalid format or contains non-numeric "
                    "characters."
                );
    }

    return decVal;
}

size_t MeshReader::safeFluentIndexConvert
(
    size_t fluentIdx,
    const std::string& context
)
{
    if (fluentIdx == 0)
    {
        throw   std::runtime_error
                (
                    "Error: Invalid Fluent index 0 in " + context
                  + ". Fluent uses 1-based indexing."
                );
    }
    return fluentIdx - 1;
}


// *************************** Main parsing method ****************************

void MeshReader::parseFile(const std::string& filePath)
{
    nodes_.clear();
    faces_.clear();
    cells_.clear();
    boundaryPatches_.clear();

    std::string token;

    std::ifstream ifs(filePath);
    if (!ifs.is_open())
    {
        throw   std::runtime_error
                (
                    "Error: Could not open mesh file: " + filePath
                );
    }

    while (ifs >> token)
    {
        if (token == MSH_COMMENT)
        {
            parseCommentSection(ifs);
        }

        if (token == MSH_DIMENSION)
        {
            parseDimensionSection(ifs);
        }

        if (token == MSH_NODES)
        {
            ifs >> token;
            parseNodesSection(ifs, token);
        }

        if (token == MSH_CELLS)
        {
            ifs >> token;
            parseCellsSection(ifs, token);
        }

        if (token == MSH_FACES)
        {
            ifs >> token;
            parseFacesSection(ifs, token);
        }

        if (token == MSH_BOUNDARIES)
        {
            ifs >> token;
            parseBoundariesSection(ifs, token);
        }
    }

    ifs.close();

    buildTopology();
    validateMesh();
    printSummary();
}


// ************************* Section parsing methods **************************

void MeshReader::parseCommentSection(std::ifstream& ifs) const
{
    int paren_level = 1;
    char ch;

    while (ifs.get(ch))
    {
        if (ch == '(')
            paren_level++;
        else if (ch == ')')
            paren_level--;
        if (paren_level == 0)
            break;
    }
}

void MeshReader::parseDimensionSection(std::ifstream& ifs) const
{
    std::string dimension;
    ifs >> dimension;
    dimension.pop_back();
    if (dimension == "2")
    {
        throw   std::runtime_error
                (
                    "This code doesn't handle 2D geometries!"
                );
    }
}

void MeshReader::parseNodesSection
(
    std::ifstream& ifs,
    const std::string& token
)
{
    if (token == "(0")
    {
        std::string firstIdxStr, lastIdxStr, zeroStrOrType;
        ifs >> firstIdxStr;
        ifs >> lastIdxStr;
        ifs >> zeroStrOrType;

        size_t lastIdx = hexToDec(lastIdxStr);

        nodes_.resize(lastIdx);
    }
    else
    {
        std::string zoneIdxStr, startIdxStr, endIdxStr;
        std::string typeStr, dimensionStr;

        // Remove leading '('
        std::string zoneToken = token;
        zoneToken.erase(0, 1);
        zoneIdxStr = zoneToken;

        ifs >> startIdxStr;
        ifs >> endIdxStr;
        ifs >> typeStr;
        ifs >> dimensionStr;

        dimensionStr.pop_back();
        dimensionStr.pop_back();

        size_t startIdx = hexToDec(startIdxStr);
        size_t endIdx = hexToDec(endIdxStr);
        size_t dimension = hexToDec(dimensionStr);

        size_t nodeGlobalIdx =
            safeFluentIndexConvert(startIdx, "node start index");

        for (size_t i = startIdx; i <= endIdx; ++i)
        {
            if (nodeGlobalIdx < endIdx && nodeGlobalIdx < nodes_.size())
            {
                Scalar xVal, yVal, zVal = 0.0;

                ifs >> xVal;
                ifs >> yVal;

                if (dimension == 3)
                {
                    ifs >> zVal;
                }

                nodes_[nodeGlobalIdx].setX(xVal);
                nodes_[nodeGlobalIdx].setY(yVal);
                nodes_[nodeGlobalIdx].setZ(zVal);

                nodeGlobalIdx++;
            }
            else
            {
                throw   std::runtime_error
                        (
                            "Error: Node index " 
                          + std::to_string(nodeGlobalIdx)
                          + " exceeds allocated node vector size "
                          + std::to_string(nodes_.size())
                        );
            }
        }
    }
}

void MeshReader::parseCellsSection
(
    std::ifstream& ifs,
    const std::string& token
)
{
    if (token == "(0")
    {
        std::string firstIdxStr, lastIdxStr, zeroStrOrType;
        ifs >> firstIdxStr;
        ifs >> lastIdxStr;
        ifs >> zeroStrOrType;

        size_t lastIdx = hexToDec(lastIdxStr);

        cells_.resize(lastIdx);
    }
}

void MeshReader::parseFacesSection
(
    std::ifstream& ifs,
    const std::string& token
)
{
    if (token == "(0")
    {
        std::string firstIdxStr, lastIdxStr, zeroStrOrType;
        ifs >> firstIdxStr;
        ifs >> lastIdxStr;
        ifs >> zeroStrOrType;

        size_t lastIdx = hexToDec(lastIdxStr);

        faces_.resize(lastIdx);
    }
    else
    {
        std::string zoneIdxStr, startIdxStr, endIdxStr;
        std::string typeStr, elementTypeStr;

        // Remove leading '('
        std::string zoneToken = token;
        zoneToken.erase(0, 1);
        zoneIdxStr = zoneToken;

        ifs >> startIdxStr;
        ifs >> endIdxStr;
        ifs >> typeStr;
        ifs >> elementTypeStr;

        elementTypeStr.pop_back();
        elementTypeStr.pop_back();

        // Check mixed element section (node count included)
        bool hasMixedElements = (elementTypeStr == "0");

        size_t zoneIdx = hexToDec(zoneIdxStr);
        size_t startIdx = hexToDec(startIdxStr);
        size_t endIdx = hexToDec(endIdxStr);

        boundaryPatches_.emplace_back(
            zoneIdx,
            safeFluentIndexConvert(startIdx, "face zone start index"),
            safeFluentIndexConvert(endIdx, "face zone end index")
        );

        std::string line;
        std::getline(ifs, line);

        for (size_t faceIdx = startIdx; faceIdx <= endIdx; ++faceIdx)
        {
            if ((faceIdx - 1) >= faces_.size())
            {
                throw   std::runtime_error
                        (
                            "Error: Face index " 
                          + std::to_string((faceIdx - 1))
                          + " is out of bounds for faces vector of "
                          + std::to_string(faces_.size()) + "."
                        );
            }

            Face &currentFace = faces_[(faceIdx - 1)];
            currentFace.setIdx(faceIdx - 1);
            currentFace.clearNodeIndices();

            if (!std::getline(ifs, line))
            {
                throw   std::runtime_error
                        (
                            "Error: Could not read face data line for id: "
                          + std::to_string(faceIdx)
                          + ". End of file reached unexpectedly."
                        );
            }

            if (line.empty())
            {
                throw   std::runtime_error
                        (
                            "Error: Empty line: face id: "
                          + std::to_string(faceIdx)
                        );
            }

            std::stringstream lineStream(line);
            std::string itemHex;
            std::vector<std::string> hexItems;

            // Put hex items in a vector
            while (lineStream >> itemHex)
            {
                hexItems.push_back(itemHex);
            }

            // If mixed elements, discard the first item (node count)
            if (hasMixedElements)
            {
                hexItems.erase(hexItems.begin());
            }

            // The last two ids are for owner and neighbor cells.
            std::string neighborHex = hexItems.back();
            hexItems.pop_back();
            std::string ownerHex = hexItems.back();
            hexItems.pop_back();

            // The remaining items are node indices
            for (const std::string &nodeIdxHex : hexItems)
            {
                currentFace.addNodeIndex
                (
                    safeFluentIndexConvert
                    (
                        hexToDec(nodeIdxHex), 
                        "face node index"
                    )
                );
            }

            currentFace.setOwnerCell(
                safeFluentIndexConvert(hexToDec(ownerHex), "owner cell index")
            );

            if (neighborHex != "0")
            {
                currentFace.setNeighborCell(
                    safeFluentIndexConvert
                    (
                        hexToDec(neighborHex), 
                        "neighbor cell index"
                    )
                );
            }
            else
            {
                currentFace.setNeighborCell(std::nullopt);
            }
        }
    }
}

void MeshReader::parseBoundariesSection
(
    std::ifstream& ifs,
    const std::string& token
)
{
    // Remove leading '('
    std::string zoneToken = token;
    zoneToken.erase(0, 1);
    size_t zoneIdx = strToDec(zoneToken);

    std::string typeString;
    ifs >> typeString;

    std::string nameString;
    ifs >> nameString;
    nameString.pop_back();
    nameString.pop_back();
    nameString.pop_back();
    nameString.pop_back();

    for (size_t i = 0; i < boundaryPatches_.size(); i++)
    {
        if (zoneIdx == boundaryPatches_[i].zoneIdx())
        {
            boundaryPatches_[i].setPatchName(nameString);

            boundaryPatches_[i].setFluentType(typeString);

            BoundaryConditionType mappedType =
                mapFluentBCToEnum(typeString);

            boundaryPatches_[i].setType(mappedType);
        }
    }
}


// ************************* Post-processing methods **************************

void MeshReader::buildTopology()
{
    for (size_t i = 0; i < cells_.size(); ++i)
    {
        cells_[i].setIdx(i);
        cells_[i].clearFaceIndices();
        cells_[i].clearNeighborCellIndices();
    }

    std::vector<std::vector<size_t>> tempCellNeighbors(cells_.size());

    for (size_t faceIdx = 0; faceIdx < faces_.size(); ++faceIdx)
    {
        const Face &currentFace = faces_[faceIdx];

        if (currentFace.ownerCell() < cells_.size())
        {
            cells_[currentFace.ownerCell()].addFaceIndex(faceIdx);
            cells_[currentFace.ownerCell()].addFaceSign(1);

            // If this face has a valid neighborCell, then that neighborCell
            if
            (
                currentFace.neighborCell().has_value()
             && currentFace.neighborCell().value() < cells_.size()
            )
            {
                tempCellNeighbors[currentFace.ownerCell()].push_back
                (
                    currentFace.neighborCell().value()
                );
            }
        }

        if
        (
            currentFace.neighborCell().has_value()
         && currentFace.neighborCell().value() < cells_.size()
        )
        {
            size_t neighborIdx = currentFace.neighborCell().value();
            cells_[neighborIdx].addFaceIndex(faceIdx);
            cells_[neighborIdx].addFaceSign(-1);

            if (currentFace.ownerCell() < cells_.size())
            {
                tempCellNeighbors[neighborIdx].push_back(
                    currentFace.ownerCell());
            }
        }
    }

    for (size_t i = 0; i < cells_.size(); ++i)
    {
        if (i < cells_.size() && i < tempCellNeighbors.size())
        {
            std::vector<size_t> &neighbors = tempCellNeighbors[i];

            std::sort(neighbors.begin(), neighbors.end());

            neighbors.erase
            (
                std::unique
                (
                    neighbors.begin(),
                    neighbors.end()
                ),
                neighbors.end()
            );

            cells_[i].setNeighborCellIndices(neighbors);
        }
        else
        {
            if (i < cells_.size())
            {
                cells_[i].clearNeighborCellIndices();
            }
            std::cerr
                << "Warning: Cell with index " << i
                << " could not have its neighbors finalized due to "
                << "indexing mismatch. cells_.size(): "
                << cells_.size()
                << ", tempCellNeighbors.size(): "
                << tempCellNeighbors.size() << std::endl;
        }
    }
}

void MeshReader::validateMesh() const
{
    for (size_t i = 0; i < faces_.size(); ++i)
    {
        if (faces_[i].nodeIndices().size() < 3)
        {
            throw   std::runtime_error
                    (
                        "Error: Face "
                      + std::to_string(i)
                      + " has only "
                      + std::to_string(faces_[i].nodeIndices().size())
                      + " nodes, minimum required is 3."
                    );
        }
    }
}

void MeshReader::printSummary() const
{
    std::cout
        << "Mesh loaded successfully:" << std::endl;

    std::cout
        << "  - Nodes: " << nodes_.size() << std::endl;

    std::cout
        << "  - Faces: " << faces_.size() << std::endl;

    std::cout
        << "  - Cells: " << cells_.size() << std::endl;

    std::cout
        << "  - Boundary patches: "
        << boundaryPatches_.size() << std::endl;
}
