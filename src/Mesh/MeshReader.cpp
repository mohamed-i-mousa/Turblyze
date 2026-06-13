/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshReader.cpp
 * @brief Implementation of the MeshReader class for Fluent mesh files
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "MeshReader.h"

// Standard library headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <charconv>
#include <vector>

// Project headers
#include "ErrorHandler.h"

// ************************* Special Member Functions *************************

MeshReader::MeshReader(const FilePath& filePath)
{
    parseFile(filePath);
}

// ****************************** Private Methods *****************************

void MeshReader::parseFile(const FilePath& filePath)
{
    Token token;

    std::ifstream ifs(filePath);
    if (!ifs.is_open())
    {
        FatalError("Could not open mesh file: " + filePath);
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
    Token dimension;
    if (!(ifs >> dimension) || dimension.empty() || dimension.back() != ')')
    {
        FatalError("Malformed dimension section in mesh file.");
    }

    // Remove trailing parenthesis
    dimension.pop_back();

    if (dimension != "3")
    {
        FatalError("This code handles 3D meshes only!");
    }
}


void MeshReader::parseNodesSection
(
    std::ifstream& ifs,
    TokenRef token
)
{
    // Header with "(0"
    if (token == "(0")
    {
        Token firstIdxStr, lastIdxStr, zeroStrOrType;
        ifs >> firstIdxStr;
        ifs >> lastIdxStr;
        ifs >> zeroStrOrType;

        const Count lastIdx = hexToDec(lastIdxStr);

        nodes_.resize(lastIdx);
    }
    // Data block with "(zoneIdx"
    else
    {
        Token zoneIdxStr;
        Token startIdxStr;
        Token endIdxStr;
        Token typeStr;
        Token dimensionStr;

        Token zoneToken = Token{token};

        if (!zoneToken.starts_with('('))
        {
            FatalError
            (
                "Malformed nodes-section header: expected zone token starting "
                "with '(' but found '" + zoneToken + "'."
            );
        }

        // Remove leading '('
        zoneToken.erase(0, 1);

        zoneIdxStr = zoneToken;

        ifs >> startIdxStr;
        ifs >> endIdxStr;
        ifs >> typeStr;
        ifs >> dimensionStr;

        if (!ifs)
        {
            FatalError
            (
                "Malformed nodes-section header: unexpected end of file "
                "while reading the zone descriptor."
            );
        }

        // The dimension token carries two trailing characters ")("
        if (dimensionStr.size() < 2)
        {
            FatalError
            (
                "Malformed nodes-section header: dimension token '"
              + dimensionStr + "' is too short."
            );
        }

        if (!dimensionStr.ends_with(")("))
        {
            FatalError
            (
                "Malformed nodes-section header: dimension token '"
              + dimensionStr
              + "' does not end with expected ')('."
            );
        }

        // Remove trailing parentheses
        dimensionStr.pop_back();
        dimensionStr.pop_back();

        const Count startIdx = hexToDec(startIdxStr);
        const Count endIdx = hexToDec(endIdxStr);
        const Count dimension = hexToDec(dimensionStr);

        for (Count idx = startIdx; idx <= endIdx; ++idx)
        {
            const Index globalIdx = idx - 1;

            if (globalIdx >= nodes_.size())
            {
                FatalError
                (
                    "Node index "
                  + std::to_string(globalIdx)
                  + " exceeds allocated node vector size "
                  + std::to_string(nodes_.size())
                );
            }

            Scalar xVal = S(0.0);
            Scalar yVal = S(0.0);
            Scalar zVal = S(0.0);

            ifs >> xVal;
            ifs >> yVal;

            if (dimension == 3)
            {
                ifs >> zVal;
            }

            if (!ifs)
            {
                FatalError
                (
                    "Malformed nodes section: failed to read coordinates "
                    "for node " + std::to_string(globalIdx) + "."
                );
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
    TokenRef token
)
{
    // Only header with "(0" is expected for cells section, no data block
    if (token == "(0")
    {
        Token firstIdxStr;
        Token lastIdxStr;
        Token zeroStrOrType;
        ifs >> firstIdxStr;
        ifs >> lastIdxStr;
        ifs >> zeroStrOrType;

        const Count lastIdx = hexToDec(lastIdxStr);

        cells_.resize(lastIdx);
    }
}


void MeshReader::parseFacesSection
(
    std::ifstream& ifs,
    TokenRef token
)
{
    // The faces section header
    if (token == "(0")
    {
        Token firstIdxStr;
        Token lastIdxStr;
        Token zeroStrOrType;
        ifs >> firstIdxStr;
        ifs >> lastIdxStr;
        ifs >> zeroStrOrType;

        const Count lastIdx = hexToDec(lastIdxStr);

        faces_.resize(lastIdx);
    }
    // Data block with "(zoneIdx" for each face zone
    else
    {
        Token zoneIdxStr;
        Token startIdxStr;
        Token endIdxStr;
        Token typeStr;
        Token elementTypeStr;

        Token zoneToken = Token{token};

        // Remove leading '('
        if (!zoneToken.starts_with('('))
        {
            FatalError
            (
                "Malformed faces-section header: expected zone token starting "
                "with '(' but found '" + zoneToken + "'."
            );
        }

        zoneToken.erase(0, 1);

        zoneIdxStr = zoneToken;

        ifs >> startIdxStr;
        ifs >> endIdxStr;
        ifs >> typeStr;
        ifs >> elementTypeStr;

        // Remove trailing parentheses
        if (!elementTypeStr.ends_with(")("))
        {
            FatalError
            (
                "Malformed faces-section header: element type token '"
              + elementTypeStr
              + "' does not end with expected ')('."
            );
        }

        elementTypeStr.pop_back();
        elementTypeStr.pop_back();

        // Check mixed element section (node count included)
        const bool hasMixedElements = (elementTypeStr == "0");

        const Index zoneIdx = hexToDec(zoneIdxStr);
        const Count startIdx = hexToDec(startIdxStr);
        const Count endIdx = hexToDec(endIdxStr);

        // If not internal face zone, create a boundary patch for this zone
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

        Token line;
        std::getline(ifs, line);
        Token itemHex;
        TokenList hexItems;

        for
        (
            Count faceIdx = startIdx;
            faceIdx <= endIdx;
            ++faceIdx
        )
        {
            if ((faceIdx - 1) >= faces_.size())
            {
                FatalError
                (
                    "Face index "
                  + std::to_string(faceIdx - 1)
                  + " is out of bounds for faces vector of "
                  + std::to_string(faces_.size())
                );
            }

            const Index internalFaceIdx = faceIdx - 1;

            Face& currentFace = faces_[internalFaceIdx];
            currentFace.setIdx(internalFaceIdx);
            currentFace.clearNodeIndices();

            if (!std::getline(ifs, line))
            {
                FatalError
                (
                    "Could not read face data line for idx: "
                  + std::to_string(faceIdx)
                  + ". End of file reached unexpectedly."
                );
            }

            if (line.empty())
            {
                FatalError
                (
                    "Empty line: face idx: "
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
            const Index nodeStart = hasMixedElements ? 1 : 0;

            if (hexItems.size() < 2)
            {
                FatalError
                (
                    "Malformed face data line for face idx: "
                  + std::to_string(faceIdx)
                );
            }

            const Index nodeEnd = hexItems.size() - 2;
            const TokenRef ownerHex = hexItems[hexItems.size() - 2];
            const TokenRef neighborHex = hexItems[hexItems.size() - 1];

            for (Index j = nodeStart; j < nodeEnd; ++j)
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
    TokenRef token
)
{
    Token zoneToken = Token{token};

    // Remove leading '('
    if (!zoneToken.starts_with('('))
    {
        FatalError
        (
            "Malformed boundaries-section header: expected zone token starting "
            "with '(' but found '" + zoneToken + "'."
        );
    }

    zoneToken.erase(0, 1);

    const Index zoneIdx = strToDec(zoneToken);

    Token bcType;
    ifs >> bcType;

    Name patchName;
    ifs >> patchName;
    
    // Remove trailing parentheses
    if (!patchName.ends_with(")())"))
    {
        FatalError
        (
            "Malformed boundaries-section header: expected name token ending "
            "with ')())' but found '" + patchName + "'."
        );
    }

    patchName.pop_back();
    patchName.pop_back();
    patchName.pop_back();
    patchName.pop_back();

    for (Index i = 0; i < boundaryPatches_.size(); ++i)
    {
        if (zoneIdx == boundaryPatches_[i].zoneIdx())
        {
            boundaryPatches_[i].setPatchName(patchName);

            boundaryPatches_[i].setType(mapFluentBCToEnum(bcType));
        }
    }
}

void MeshReader::buildTopology()
{
    for (Index cellIdx = 0; cellIdx < cells_.size(); ++cellIdx)
    {
        cells_[cellIdx].setIdx(cellIdx);
    }

    std::vector<IndexList> tempCellNeighbors(cells_.size());

    for (Index faceIdx = 0; faceIdx < faces_.size(); ++faceIdx)
    {
        const Face& currentFace = faces_[faceIdx];

        // Map face to owner cell and assign sign +1
        if (currentFace.ownerCell() < cells_.size())
        {
            cells_[currentFace.ownerCell()].addFace(faceIdx, 1);

            // If this face has a valid neighborCell
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

        // If this face has a valid neighborCell, 
        // map face to neighbor cell and assign sign -1
        if
        (
            currentFace.neighborCell().has_value()
         && currentFace.neighborCell().value() < cells_.size()
        )
        {
            const Index neighborIdx = currentFace.neighborCell().value();
            
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

    // After processing all faces, set the neighbor cell indices for each cell
    for (Index cellIdx = 0; cellIdx < cells_.size(); ++cellIdx)
    {
        cells_[cellIdx].setNeighborCellIndices(tempCellNeighbors[cellIdx]);
    }
}


void MeshReader::validateMesh() const
{
    for (Index faceIdx = 0; faceIdx < faces_.size(); ++faceIdx)
    {
        const Face& face = faces_[faceIdx];

        if (face.nodeIndices().size() < 3)
        {
            FatalError
            (
                "Face "
              + std::to_string(faceIdx)
              + " has only "
              + std::to_string(face.nodeIndices().size())
              + " nodes, minimum required is 3."
            );
        }

        if (face.ownerCell() >= cells_.size())
        {
            FatalError
            (
                "Face "
              + std::to_string(faceIdx)
              + " has owner cell index "
              + std::to_string(face.ownerCell())
              + " outside the valid range [0, "
              + std::to_string(cells_.size())
              + ")."
            );
        }

        if
        (
            face.neighborCell().has_value()
         && face.neighborCell().value() >= cells_.size()
        )
        {
            FatalError
            (
                "Face "
              + std::to_string(faceIdx)
              + " has neighbor cell index "
              + std::to_string(face.neighborCell().value())
              + " outside the valid range [0, "
              + std::to_string(cells_.size())
              + ")."
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

PatchType MeshReader::mapFluentBCToEnum(TokenRef fluentType)
{
    for (const auto& [name, type] : bcMappings_)
    {
        if (fluentType == name) return type;
    }

    Warning
    (
        "Unknown Fluent boundary type encountered: "
      + Token(fluentType)
    );

    return PatchType::undefined;
}


Count MeshReader::hexToDec(TokenRef hexStr)
{
    Count decVal = 0;

    const auto result =
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
        FatalError
        (
            "Failed to convert hex string '"
          + Token(hexStr) + "' to decimal."
        );
    }

    return decVal;
}


Count MeshReader::strToDec(TokenRef decStr)
{
    Count decVal = 0;
    
    const auto result =
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
        FatalError
        (
            "Failed to convert decimal string '"
          + Token(decStr)
          + "' to Count. Invalid format or "
            "contains non-numeric characters."
        );
    }

    return decVal;
}


Index MeshReader::safeFluentIndexConvert
(
    Count fluentIdx,
    MessageRef context
)
{
    if (fluentIdx == 0)
    {
        FatalError
        (
            "Invalid Fluent index 0 in "
          + Message(context)
          + ". Fluent uses 1-based indexing."
        );
    }
    return fluentIdx - 1;
}
