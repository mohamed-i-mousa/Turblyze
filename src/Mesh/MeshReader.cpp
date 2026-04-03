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


// ******************************* Constructor ********************************

MeshReader::MeshReader(const std::string& filePath)
{
    parseFile(filePath);
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
        throw
            std::runtime_error
            (
                "Error: Could not open mesh file: "
              + filePath
            );
    }

    while (ifs >> token)
    {
        if (token == MSH_COMMENT)
        {
            parseCommentSection(ifs);
        }

        else if (token == MSH_DIMENSION)
        {
            parseDimensionSection(ifs);
        }

        else if (token == MSH_NODES)
        {
            ifs >> token;
            parseNodesSection(ifs, token);
        }

        else if (token == MSH_CELLS)
        {
            ifs >> token;
            parseCellsSection(ifs, token);
        }

        else if (token == MSH_FACES)
        {
            ifs >> token;
            parseFacesSection(ifs, token);
        }

        else if (token == MSH_BOUNDARIES)
        {
            ifs >> token;
            parseBoundariesSection(ifs, token);
        }
    }

    buildTopology();
    validateMesh();
    printSummary();
}


// ************************* Section parsing methods **************************

void MeshReader::parseCommentSection(std::ifstream& ifs) const
{
    int parenLevel = 1;
    bool inQuotes = false;
    char ch;

    while (ifs.get(ch))
    {
        if (ch == '"')
            inQuotes = !inQuotes;

        if (!inQuotes)
        {
            if (ch == '(')
                parenLevel++;
            else if (ch == ')')
                parenLevel--;
        }

        if (parenLevel == 0)
            break;
    }
}

void MeshReader::parseDimensionSection(std::ifstream& ifs) const
{
    std::string dimension;
    ifs >> dimension;
    dimension.pop_back();
    if (dimension != "3")
    {
        throw
            std::runtime_error
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

        for (size_t idx = startIdx; idx <= endIdx; ++idx)
        {
            size_t globalIdx = idx - 1;

            if (globalIdx >= nodes_.size())
            {
                throw
                    std::runtime_error
                    (
                        "Error: Node index "
                      + std::to_string(globalIdx)
                      + " exceeds allocated node vector"
                      + " size "
                      + std::to_string(nodes_.size())
                    );
            }

            Scalar xVal = 0.0;
            Scalar yVal = 0.0;
            Scalar zVal = 0.0;

            ifs >> xVal;
            ifs >> yVal;

            if (dimension == 3)
            {
                ifs >> zVal;
            }

            nodes_[globalIdx].setX(xVal);
            nodes_[globalIdx].setY(yVal);
            nodes_[globalIdx].setZ(zVal);
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

        if (typeStr != "2")
        {
            boundaryPatches_.emplace_back
            (
                zoneIdx,
                safeFluentIndexConvert
                (
                    startIdx, "face zone start index"
                ),
                safeFluentIndexConvert
                (
                    endIdx, "face zone end index"
                )
            );
        }

        std::string line;
        std::getline(ifs, line);

        std::string itemHex;
        std::vector<std::string> hexItems;

        for
        (
            size_t faceIdx = startIdx;
            faceIdx <= endIdx;
            ++faceIdx
        )
        {
            if ((faceIdx - 1) >= faces_.size())
            {
                throw
                    std::runtime_error
                    (
                        "Error: Face index "
                      + std::to_string((faceIdx - 1))
                      + " is out of bounds for faces"
                      + " vector of "
                      + std::to_string(faces_.size())
                      + "."
                    );
            }

            Face &currentFace = faces_[(faceIdx - 1)];
            currentFace.setIdx(faceIdx - 1);
            currentFace.clearNodeIndices();

            if (!std::getline(ifs, line))
            {
                throw
                    std::runtime_error
                    (
                        "Error: Could not read face"
                      + std::string(" data line for id: ")
                      + std::to_string(faceIdx)
                      + ". End of file reached"
                      + " unexpectedly."
                    );
            }

            if (line.empty())
            {
                throw
                    std::runtime_error
                    (
                        "Error: Empty line: face id: "
                      + std::to_string(faceIdx)
                    );
            }

            std::istringstream lineStream(line);
            hexItems.clear();

            // Put hex items in a vector
            while (lineStream >> itemHex)
            {
                hexItems.push_back(itemHex);
            }

            // Owner and neighbor are the last two items;
            // node indices are everything before them.
            // If mixed elements, the first item is the
            // node count and must be skipped.
            size_t nodeStart = hasMixedElements ? 1 : 0;
            size_t nodeEnd   = hexItems.size() - 2;

            std::string_view ownerHex =
                hexItems[hexItems.size() - 2];
            std::string_view neighborHex =
                hexItems[hexItems.size() - 1];

            for (size_t j = nodeStart; j < nodeEnd; ++j)
            {
                currentFace.addNodeIndex
                (
                    safeFluentIndexConvert
                    (
                        hexToDec(hexItems[j]),
                        "face node index"
                    )
                );
            }

            currentFace.setOwnerCell
            (
                safeFluentIndexConvert
                (
                    hexToDec(ownerHex),
                    "owner cell index"
                )
            );

            if (neighborHex != "0")
            {
                currentFace.setNeighborCell
                (
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

            boundaryPatches_[i].setType
            (
                mapFluentBCToEnum(typeString)
            );
        }
    }
}


// ************************* Post-processing methods **************************

void MeshReader::buildTopology()
{
    for (size_t i = 0; i < cells_.size(); ++i)
    {
        cells_[i].setIdx(i);
    }

    std::vector<std::vector<size_t>> tempCellNeighbors(cells_.size());

    for (size_t faceIdx = 0; faceIdx < faces_.size(); ++faceIdx)
    {
        const Face &currentFace = faces_[faceIdx];

        if (currentFace.ownerCell() < cells_.size())
        {
            cells_[currentFace.ownerCell()].addFace(faceIdx, 1);

            // If this face has a valid neighborCell
            if
            (
                currentFace.neighborCell().has_value()
             && currentFace.neighborCell().value()
              < cells_.size()
            )
            {
                tempCellNeighbors
                    [currentFace.ownerCell()].push_back
                (
                    currentFace.neighborCell().value()
                );
            }
        }

        if
        (
            currentFace.neighborCell().has_value()
         && currentFace.neighborCell().value()
          < cells_.size()
        )
        {
            size_t neighborIdx = currentFace.neighborCell().value();
            
            cells_[neighborIdx].addFace(faceIdx, -1);

            if (currentFace.ownerCell() < cells_.size())
            {
                tempCellNeighbors[neighborIdx].push_back
                (
                    currentFace.ownerCell()
                );
            }
        }
    }

    for (size_t i = 0; i < cells_.size(); ++i)
    {
        std::vector<size_t> &neighbors =
            tempCellNeighbors[i];

        std::sort
        (
            neighbors.begin(), neighbors.end()
        );

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
}

void MeshReader::validateMesh() const
{
    for (size_t i = 0; i < faces_.size(); ++i)
    {
        if (faces_[i].nodeIndices().size() < 3)
        {
            throw
                std::runtime_error
                (
                    "Error: Face "
                  + std::to_string(i)
                  + " has only "
                  + std::to_string
                    (
                        faces_[i].nodeIndices().size()
                    )
                  + " nodes, minimum required is 3."
                );
        }
    }
}

void MeshReader::printSummary() const
{
    std::cout
        << "Mesh loaded successfully:" << '\n';

    std::cout
        << "  - Nodes: " << nodes_.size() << '\n';

    std::cout
        << "  - Faces: " << faces_.size() << '\n';

    std::cout
        << "  - Cells: " << cells_.size() << '\n';

    std::cout
        << "  - Boundary patches: "
        << boundaryPatches_.size() << '\n';
}


// ************************* Static utility methods ***************************

PatchType MeshReader::mapFluentBCToEnum
(
    std::string_view fluentType
)
{
    for (const auto& [name, type] : bcMappings_)
    {
        if (fluentType == name)
            return type;
    }

    std::cerr
        << "Warning: Unknown Fluent boundary type encountered: "
        << fluentType << std::endl;

    return PatchType::UNDEFINED;
}

size_t MeshReader::hexToDec(std::string_view hexStr)
{
    size_t decVal;
    auto result =
        std::from_chars
        (
            hexStr.data(),
            hexStr.data() + hexStr.size(),
            decVal,
            16
        );

    if
    (
        result.ec != std::errc()
     || result.ptr != hexStr.data() + hexStr.size()
    )
    {
        throw
            std::runtime_error
            (
                "Error: Failed to convert hex string"
              + std::string(" '") + std::string(hexStr)
              + "' to decimal."
            );
    }

    return decVal;
}

size_t MeshReader::strToDec(std::string_view decStr)
{
    if (decStr.empty())
    {
        throw
            std::runtime_error
            (
                "Error: Attempted to convert empty"
              + std::string(" decimal string to ")
              + "size_t."
            );
    }

    size_t decVal;
    auto result =
        std::from_chars
        (
            decStr.data(),
            decStr.data() + decStr.size(),
            decVal
        );

    if
    (
        result.ec != std::errc()
     || result.ptr != decStr.data() + decStr.size()
    )
    {
        throw
            std::runtime_error
            (
                "Error: Failed to convert decimal"
              + std::string(" string '") + std::string(decStr)
              + "' to size_t. Invalid format or"
              + " contains non-numeric characters."
            );
    }

    return decVal;
}

size_t MeshReader::safeFluentIndexConvert
(
    size_t fluentIdx,
    std::string_view context
)
{
    if (fluentIdx == 0)
    {
        throw
            std::runtime_error
            (
                "Error: Invalid Fluent index 0 in "
              + std::string(context)
              + ". Fluent uses 1-based indexing."
            );
    }
    return fluentIdx - 1;
}
