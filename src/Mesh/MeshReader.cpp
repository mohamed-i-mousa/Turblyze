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
#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <format>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string_view>
#include <vector>

// Project headers
#include "ErrorHandler.h"

// ***************************** Internal Helpers *****************************

namespace
{

/// Fluent mesh-file section identifier tokens (ASCII data blocks)
constexpr std::string_view MSH_COMMENT    = "(0";
constexpr std::string_view MSH_DIMENSION  = "(2";
constexpr std::string_view MSH_NODES      = "(10";
constexpr std::string_view MSH_CELLS      = "(12";
constexpr std::string_view MSH_FACES      = "(13";
constexpr std::string_view MSH_BOUNDARIES = "(45";

/// Fluent binary data-block identifiers (Binary data blocks)
constexpr std::string_view MSH_NODES_SP   = "(2010";
constexpr std::string_view MSH_NODES_DP   = "(3010";
constexpr std::string_view MSH_CELLS_SP   = "(2012";
constexpr std::string_view MSH_CELLS_DP   = "(3012";
constexpr std::string_view MSH_FACES_SP   = "(2013";
constexpr std::string_view MSH_FACES_DP   = "(3013";

// Convert hexadecimal string to Count
[[nodiscard]] Count hexToDec(const Token& hexStr)
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


// Convert decimal string to Count
[[nodiscard]] Count strToDec(const Token& decStr)
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


// Safely convert 1-based Fluent index to 0-based index
[[nodiscard]] Index safeFluentIndexConvert
(
    Count fluentIdx,
    Message context
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


// Strip the single trailing ')' that closes a binary section header
void stripBinaryHeaderClose(Token& token)
{
    if (token.empty() || token.back() != ')')
    {
        FatalError
        (
            "Malformed binary section header: token '"
          + token + "' does not end with expected ')'."
        );
    }
    token.pop_back();
}


// Advance the stream past the '(' that introduces a raw data block
void advanceToRawData(std::ifstream& ifs)
{
    char ch = '\0';
    while (ifs.get(ch))
    {
        if (ch == '(')
        {
            return;
        }
    }

    FatalError("Malformed binary block: missing '(' before raw data.");
}


// Read an exact number of raw bytes into buffer, guarding against a short read
void readRawBytes
(
    std::ifstream& ifs,
    char* buffer,
    Count byteCount
)
{
    ifs.read(buffer, static_cast<std::streamsize>(byteCount));

    if (!ifs || static_cast<Count>(ifs.gcount()) != byteCount)
    {
        FatalError
        (
            "Unexpected end of file while reading a binary data block."
        );
    }
}


// Read one little-endian 32-bit integer from the raw data block
[[nodiscard]] std::int32_t readInt32(std::ifstream& ifs)
{
    std::int32_t value = 0;
    ifs.read(reinterpret_cast<char*>(&value), sizeof(std::int32_t));

    if (!ifs)
    {
        FatalError
        (
            "Unexpected end of file while reading binary connectivity."
        );
    }

    return value;
}


// Consume the "\n)\nEnd of Binary Section <id>)" trailer
void consumeBinaryTrailer(std::ifstream& ifs)
{
    Token token;
    while (ifs >> token)
    {
        if (token == "Section")
        {
            // The following token is the "<id>)" closing marker
            ifs >> token;
            return;
        }
    }

    FatalError
    (
        "Malformed binary block: reached end of file while scanning for "
        "the 'End of Binary Section' trailer."
    );
}

} // namespace

// ************************* Special Member Functions *************************

MeshReader::MeshReader(const FilePath& filePath)
{
    parseFile(filePath);
}

// ****************************** Private Methods *****************************

void MeshReader::parseFile(const FilePath& filePath)
{
    Token token;

    std::ifstream ifs(filePath, std::ios::binary);
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

        else if (token == MSH_NODES_DP || token == MSH_NODES_SP)
        {
            const bool doublePrecision = (token == MSH_NODES_DP);
            ifs >> token;
            parseNodesSectionBinary(ifs, token, doublePrecision);
        }

        else if (token == MSH_CELLS)
        {
            ifs >> token;
            parseCellsSection(ifs, token);
        }

        else if (token == MSH_CELLS_DP || token == MSH_CELLS_SP)
        {
            ifs >> token;
            parseCellsSectionBinary(ifs, token);
        }

        else if (token == MSH_FACES)
        {
            ifs >> token;
            parseFacesSection(ifs, token);
        }

        else if (token == MSH_FACES_DP || token == MSH_FACES_SP)
        {
            ifs >> token;
            parseFacesSectionBinary(ifs, token);
        }

        else if (token == MSH_BOUNDARIES)
        {
            ifs >> token;
            parseBoundariesSection(ifs, token);
        }
    }

    buildTopology();
    validateMesh();

    std::cout << std::format
    (
        "Mesh loaded successfully:\n"
        "  - Nodes: {}\n"
        "  - Faces: {}\n"
        "  - Cells: {}\n"
        "  - Boundary patches: {}\n",
        nodes_.size(), faces_.size(), cells_.size(), boundaryPatches_.size()
    );
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
    const Token& token
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

            nodes_[globalIdx] = Vector(xVal, yVal, zVal);
        }
    }
}


void MeshReader::parseNodesSectionBinary
(
    std::ifstream& ifs,
    const Token& token,
    bool doublePrecision
)
{
    Token zoneToken = Token{token};

    if (!zoneToken.starts_with('('))
    {
        FatalError
        (
            "Malformed binary nodes-section header: expected zone token "
            "starting with '(' but found '" + zoneToken + "'."
        );
    }

    Token startIdxStr;
    Token endIdxStr;
    Token typeStr;
    Token dimensionStr;

    ifs >> startIdxStr;
    ifs >> endIdxStr;
    ifs >> typeStr;
    ifs >> dimensionStr;

    if (!ifs)
    {
        FatalError
        (
            "Malformed binary nodes-section header: unexpected end of file "
            "while reading the zone descriptor."
        );
    }

    stripBinaryHeaderClose(dimensionStr);

    const Index startIdx = hexToDec(startIdxStr);
    const Index endIdx = hexToDec(endIdxStr);
    const Index dimension = hexToDec(dimensionStr);

    advanceToRawData(ifs);

    const Index nNodes = endIdx - startIdx + 1;
    const Count width =
        doublePrecision ? sizeof(double) : sizeof(float);
    const Count byteCount = nNodes * dimension * width;

    std::vector<char> buffer(byteCount);
    readRawBytes(ifs, buffer.data(), byteCount);

    const char* ptr = buffer.data();

    for (Index idx = startIdx; idx <= endIdx; ++idx)
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

        Scalar coords[3] = {S(0.0), S(0.0), S(0.0)};

        for (Index component = 0; component < dimension; ++component)
        {
            if (doublePrecision)
            {
                double value = 0.0;
                std::memcpy(&value, ptr, sizeof(double));
                ptr += sizeof(double);
                coords[component] = static_cast<Scalar>(value);
            }
            else
            {
                float value = 0.0f;
                std::memcpy(&value, ptr, sizeof(float));
                ptr += sizeof(float);
                coords[component] = static_cast<Scalar>(value);
            }
        }

        nodes_[globalIdx] = Vector(coords[0], coords[1], coords[2]);
    }

    consumeBinaryTrailer(ifs);
}


void MeshReader::parseCellsSection
(
    std::ifstream& ifs,
    const Token& token
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


void MeshReader::parseCellsSectionBinary
(
    std::ifstream& ifs,
    const Token& token
)
{
    Token zoneToken = Token{token};

    if (!zoneToken.starts_with('('))
    {
        FatalError
        (
            "Malformed binary cells-section header: expected zone token "
            "starting with '(' but found '" + zoneToken + "'."
        );
    }

    Token startIdxStr;
    Token endIdxStr;
    Token typeStr;
    Token elementTypeStr;

    ifs >> startIdxStr;
    ifs >> endIdxStr;
    ifs >> typeStr;
    ifs >> elementTypeStr;

    if (!ifs)
    {
        FatalError
        (
            "Malformed binary cells-section header: unexpected end of file "
            "while reading the zone descriptor."
        );
    }

    stripBinaryHeaderClose(elementTypeStr);

    // A uniform-type cell zone carries no binary body to consume
    if (elementTypeStr != "0")
    {
        return;
    }

    const Count startIdx = hexToDec(startIdxStr);
    const Count endIdx = hexToDec(endIdxStr);
    const Count nCells = endIdx - startIdx + 1;
    const std::size_t byteCount = nCells * sizeof(std::int32_t);

    advanceToRawData(ifs);

    std::vector<char> scratch(byteCount);
    readRawBytes(ifs, scratch.data(), byteCount);

    consumeBinaryTrailer(ifs);
}


void MeshReader::parseFacesSection
(
    std::ifstream& ifs,
    const Token& token
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

        // Mixed (0) and polygonal (5) zones include a per-face node count
        const bool hasNodeCountPrefix =
            (elementTypeStr == "0" || elementTypeStr == "5");

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

            // Owner and neighbor are the last two items.
            // A leading node count, if present, must be skipped
            const Index nodeStart = hasNodeCountPrefix ? 1 : 0;

            if (hexItems.size() < 2)
            {
                FatalError
                (
                    "Malformed face data line for face idx: "
                  + std::to_string(faceIdx)
                );
            }

            const Index nodeEnd = hexItems.size() - 2;
            const Token& ownerHex = hexItems[hexItems.size() - 2];
            const Token& neighborHex = hexItems[hexItems.size() - 1];

            IndexList nodeIndices;
            if (nodeEnd > nodeStart)
            {
                nodeIndices.reserve(nodeEnd - nodeStart);
            }

            for (Index j = nodeStart; j < nodeEnd; ++j)
            {
                nodeIndices.push_back
                (
                    safeFluentIndexConvert
                    (
                        hexToDec(hexItems[j]),
                        "face node index"
                    )
                );
            }

            const Index ownerCell = safeFluentIndexConvert
            (
                hexToDec(ownerHex),
                "owner cell index"
            );

            if (neighborHex != "0")
            {
                const Index neighborCell = safeFluentIndexConvert
                (
                    hexToDec(neighborHex),
                    "neighbor cell index"
                );

                faces_[internalFaceIdx] = Face
                (
                    internalFaceIdx,
                    std::move(nodeIndices),
                    ownerCell,
                    neighborCell
                );
            }
            else
            {
                faces_[internalFaceIdx] = Face
                (
                    internalFaceIdx,
                    std::move(nodeIndices),
                    ownerCell
                );
            }
        }
    }
}


void MeshReader::parseFacesSectionBinary
(
    std::ifstream& ifs,
    const Token& token
)
{
    Token zoneToken = Token{token};

    if (!zoneToken.starts_with('('))
    {
        FatalError
        (
            "Malformed binary faces-section header: expected zone token "
            "starting with '(' but found '" + zoneToken + "'."
        );
    }

    // Remove leading '('
    zoneToken.erase(0, 1);

    Token startIdxStr;
    Token endIdxStr;
    Token typeStr;
    Token elementTypeStr;

    ifs >> startIdxStr;
    ifs >> endIdxStr;
    ifs >> typeStr;
    ifs >> elementTypeStr;

    if (!ifs)
    {
        FatalError
        (
            "Malformed binary faces-section header: unexpected end of file "
            "while reading the zone descriptor."
        );
    }

    stripBinaryHeaderClose(elementTypeStr);

    // Mixed (0) and polygonal (5) zones prefix each face record with its
    // node count; uniform zones (triangle 3, quad 4) carry a fixed count.
    const bool hasNodeCountPrefix =
        (elementTypeStr == "0" || elementTypeStr == "5");
    const Count fixedNodeCount =
        hasNodeCountPrefix ? 0 : hexToDec(elementTypeStr);

    const Index zoneIdx = hexToDec(zoneToken);
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

    advanceToRawData(ifs);

    for (Count faceIdx = startIdx; faceIdx <= endIdx; ++faceIdx)
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

        // Node count, node indices, then owner and neighbor cells
        const Count nodeCount = hasNodeCountPrefix
            ? static_cast<Count>(readInt32(ifs))
            : fixedNodeCount;

        IndexList nodeIndices;
        nodeIndices.reserve(nodeCount);

        for (Count j = 0; j < nodeCount; ++j)
        {
            nodeIndices.push_back
            (
                safeFluentIndexConvert
                (
                    static_cast<Count>(readInt32(ifs)),
                    "face node index"
                )
            );
        }

        const std::int32_t owner = readInt32(ifs);
        const std::int32_t neighbor = readInt32(ifs);

        const Index ownerCell = safeFluentIndexConvert
        (
            static_cast<Count>(owner),
            "owner cell index"
        );

        if (neighbor != 0)
        {
            const Index neighborCell = safeFluentIndexConvert
            (
                static_cast<Count>(neighbor),
                "neighbor cell index"
            );

            faces_[internalFaceIdx] = Face
            (
                internalFaceIdx,
                std::move(nodeIndices),
                ownerCell,
                neighborCell
            );
        }
        else
        {
            faces_[internalFaceIdx] = Face
            (
                internalFaceIdx,
                std::move(nodeIndices),
                ownerCell
            );
        }
    }

    consumeBinaryTrailer(ifs);
}


void MeshReader::parseBoundariesSection
(
    std::ifstream& ifs,
    const Token& token
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
            boundaryPatches_[i].setName(patchName);
            boundaryPatches_[i].setType(PatchType::physical);
        }
    }
}

void MeshReader::buildTopology()
{
    const Count numCells = cells_.size();

    std::vector<IndexList> cellFaces(numCells);
    std::vector<Cell::FaceSignList> cellSigns(numCells);
    std::vector<IndexList> cellNeighbors(numCells);

    for (Index faceIdx = 0; faceIdx < faces_.size(); ++faceIdx)
    {
        const Face& currentFace = faces_[faceIdx];

        // Map face to owner cell and assign sign +1
        if (currentFace.ownerCell() < numCells)
        {
            cellFaces[currentFace.ownerCell()].push_back(faceIdx);
            cellSigns[currentFace.ownerCell()].push_back(static_cast<std::int8_t>(1));

            // If this face has a valid neighborCell
            if
            (
                currentFace.neighborCell().has_value()
             && currentFace.neighborCell().value() < numCells
            )
            {
                cellNeighbors[currentFace.ownerCell()].push_back
                (
                    currentFace.neighborCell().value()
                );
            }
        }

        // Map face to neighbor cell and assign sign -1 if neighborCell is valid
        if
        (
            currentFace.neighborCell().has_value()
         && currentFace.neighborCell().value() < numCells
        )
        {
            const Index neighborIdx = currentFace.neighborCell().value();

            cellFaces[neighborIdx].push_back(faceIdx);
            cellSigns[neighborIdx].push_back(static_cast<std::int8_t>(-1));

            if (currentFace.ownerCell() < numCells)
            {
                cellNeighbors[neighborIdx].push_back
                (
                    currentFace.ownerCell()
                );
            }
        }
    }

    // Construct cells with their connectivity
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        cells_[cellIdx] = Cell
        (
            cellIdx,
            std::move(cellFaces[cellIdx]),
            std::move(cellNeighbors[cellIdx]),
            std::move(cellSigns[cellIdx])
        );
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
