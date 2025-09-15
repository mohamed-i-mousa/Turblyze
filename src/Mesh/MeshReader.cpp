/******************************************************************************
 * @file MeshReader.cpp
 * @brief Implementation of Fluent mesh file reader
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <charconv>

#include "MeshReader.hpp"

/// Fluent file section identifiers
const std::string MSH_COMMENT = "(0";
const std::string MSH_DIMENSION = "(2";
const std::string MSH_NODES = "(10";
const std::string MSH_CELLS = "(12";
const std::string MSH_FACES = "(13";
const std::string MSH_BOUNDARIES = "(45";

/**
 * @brief Convert hexadecimal string to size_t
 * @param hexStr Hexadecimal string to convert
 * @return Converted decimal value
 * @throws std::runtime_error if conversion fails
 */
size_t hexToDec(const std::string &hexStr) 
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
        throw std::runtime_error
            (
                "Error: Failed to convert hex string '" 
              + hexStr + "' to decimal."
            );
    }
    
    return decVal;
}

/**
 * @brief Convert decimal string to size_t
 * @param decStr Decimal string to convert
 * @return Converted decimal value
 * @throws std::runtime_error if string is empty or conversion fails
 */
size_t strToDec(const std::string& decStr) 
{
    if (decStr.empty()) 
    {
        throw std::runtime_error
        (
            "Error: Attempted to convert empty decimal string to size_t."
        );
    }

    std::stringstream ss;
    ss  << decStr; 
    size_t decVal;
    ss  >> decVal;
    
    if (ss.fail() || !ss.eof()) 
    {
        throw std::runtime_error
        (
            "Error: Failed to convert decimal string '" 
          + decStr 
          + "' to size_t. Invalid format or contains non-numeric characters."
        );
    }

    return decVal;
}

void readMshFile
(
    const std::string& filePath,
    std::vector<Vector>& allNodes,
    std::vector<Face>& allFaces,
    std::vector<Cell>& allCells,
    std::vector<BoundaryPatch>& allBoundaryPatches
)
{
    allNodes.clear();
    allFaces.clear();
    allCells.clear();
    allBoundaryPatches.clear();

    std::string line;
    std::string token;
    char ch;

    size_t numCells = 0;

    std::ifstream ifs(filePath);
    if (!ifs.is_open()) 
    {
        throw std::runtime_error
            (
                "Error: Could not open mesh file: " + filePath
            );
    }

    while (ifs >> token) 
    {
        if (token == MSH_COMMENT) 
        { 
            int paren_level = 1;

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

        if (token == MSH_DIMENSION) 
        { 
            std::string dimension;
            ifs >> dimension;
            dimension.pop_back();
            if (dimension == "2") 
            {
                throw std::runtime_error
                (
                    "This code doesn't handle 2D geometries!"
                );
            }
        }

        if (token == MSH_NODES) 
        { 
            ifs >> token;

            if (token == "(0") 
            {
                std::string firstIdxStr, lastIdxStr, zeroStrOrType;
                ifs >> firstIdxStr;
                ifs >> lastIdxStr;
                ifs >> zeroStrOrType;

                size_t lastIdx = hexToDec(lastIdxStr);

                allNodes.resize(lastIdx);
            }
            else 
            {
                std::string zoneIdxStr, startIdxStr, endIdxStr;
                std::string typeStr, dimensionStr;

                // Remove leading '('
                token.erase(0, 1);
                zoneIdxStr = token;

                ifs >> startIdxStr;
                ifs >> endIdxStr;
                ifs >> typeStr;
                ifs >> dimensionStr;

                dimensionStr.pop_back();
                dimensionStr.pop_back();

                size_t startIdx = hexToDec(startIdxStr);
                size_t endIdx = hexToDec(endIdxStr);
                size_t dimension = hexToDec(dimensionStr);

                size_t nodeGlobalIdx = startIdx - 1;

                for (size_t i = startIdx; i <= endIdx; ++i) 
                {
                    if (nodeGlobalIdx < endIdx) 
                    {
                        Scalar xVal, yVal, zVal = 0.0;

                        ifs >> xVal;
                        ifs >> yVal;

                        if (dimension == 3) 
                        {
                            ifs >> zVal;
                        }

                        allNodes[nodeGlobalIdx].setX(xVal);
                        allNodes[nodeGlobalIdx].setY(yVal);
                        allNodes[nodeGlobalIdx].setZ(zVal);
                        
                        nodeGlobalIdx++;
                    } 
                    else 
                    {
                        throw std::runtime_error
                            (
                                "Error: Node index " 
                              + std::to_string(nodeGlobalIdx)
                              + " exceeds allocated node vector size "
                              + std::to_string(allNodes.size())
                            );            
                    }
                }
            }
    }

        if (token == MSH_CELLS)
        {
            ifs >> token;
            
            if (token == "(0")
            {
                std::string firstIdxStr, lastIdxStr, zeroStrOrType;
                ifs >> firstIdxStr;
                ifs >> lastIdxStr;
                ifs >> zeroStrOrType;

                size_t lastIdx = hexToDec(lastIdxStr);

                allCells.resize(lastIdx);
                numCells = allCells.size();
            }
        }

        if (token == MSH_FACES)
        {
            ifs >> token;

            if (token == "(0")
            {
                std::string firstIdxStr, lastIdxStr, zeroStrOrType;
                ifs >> firstIdxStr;
                ifs >> lastIdxStr;
                ifs >> zeroStrOrType;

                size_t lastIdx = hexToDec(lastIdxStr);

                allFaces.resize(lastIdx);
            }
            else 
            {
                std::string zoneIdxStr, startIdxStr, endIdxStr;
                std::string typeStr, elementTypeStr;

                // Remove leading '('
                token.erase(0, 1);
                zoneIdxStr = token;

                ifs >> startIdxStr;
                ifs >> endIdxStr;
                ifs >> typeStr;
                ifs >> elementTypeStr;

                elementTypeStr.pop_back();
                elementTypeStr.pop_back();
                
                // Check if this is a mixed element section (node count included)
                bool hasMixedElements = (elementTypeStr == "0");

                size_t zoneIdx = hexToDec(zoneIdxStr);
                size_t startIdx = hexToDec(startIdxStr);
                size_t endIdx = hexToDec(endIdxStr);
                
                allBoundaryPatches.emplace_back(zoneIdx, startIdx-1, endIdx-1);

                std::string line;
                std::getline(ifs, line); 

                for (size_t faceIdx = startIdx; faceIdx <= endIdx; ++faceIdx)
                {
                    if ((faceIdx - 1) >= allFaces.size())
                    {
                        throw std::runtime_error
                            (
                                "Error: Face index " 
                              + std::to_string((faceIdx - 1))
                              + " is out of bounds for allFaces vector of "
                              + std::to_string(allFaces.size()) + "."
                            );
                    }

                    Face &currentFace = allFaces[(faceIdx - 1)];
                    currentFace.setId(faceIdx - 1);
                    currentFace.clearNodeIndices();

                    if (!std::getline(ifs, line))
                    {
                        throw std::runtime_error
                            (
                                "Error: Could not read face data line for id: "
                              + std::to_string(faceIdx)
                              + ". End of file reached unexpectedly."
                            );
                    }
                    
                    if (line.empty())
                    {
                        throw std::runtime_error
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
                            hexToDec(nodeIdxHex) - 1
                        );
                    }

                    currentFace.setOwnerCell(hexToDec(ownerHex) - 1);

                    if (neighborHex != "0")
                    {
                        currentFace.setNeighborCell(
                            hexToDec(neighborHex) - 1);
                    }
                    else 
                    {
                        currentFace.setNeighborCell(std::nullopt);
                    }
                } 
            }
    }

        if (token == MSH_BOUNDARIES)
        {
            ifs >> token;

            // Remove leading '('
            token.erase(0, 1); 
            size_t zoneIdx = strToDec(token);

            std::string typeString;
            ifs >> typeString;

            std::string nameString;
            ifs >> nameString;
            nameString.pop_back();
            nameString.pop_back();
            nameString.pop_back();
            nameString.pop_back();

            for (size_t i = 0; i < allBoundaryPatches.size(); i++)
            {
                if (zoneIdx == allBoundaryPatches[i].zoneID())
                {
                    allBoundaryPatches[i].setPatchName(nameString);
                    allBoundaryPatches[i].setFluentType(typeString);
                    BoundaryConditionType mappedType =
                        mapFluentBCToEnum(typeString);
                    allBoundaryPatches[i].setType(mappedType);
                }
            }
        }
    }

    ifs.close();

    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        allCells[i].setId(i);
        allCells[i].clearFaceIndices();
        allCells[i].clearNeighborCellIndices();
    }

    std::vector<std::vector<size_t>> tempCellNeighbors(numCells);

    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
    {
        const Face &currentFace = allFaces[faceIdx];

        if (currentFace.ownerCell() < allCells.size()) 
        {
            allCells[currentFace.ownerCell()].addFaceIndex(faceIdx);
            allCells[currentFace.ownerCell()].addFaceSign(1);

            // If this face has a valid neighborCell, then that neighborCell
            if 
            (
                currentFace.neighborCell().has_value() 
             && currentFace.neighborCell().value() < allCells.size()
            ) 
            {
                tempCellNeighbors[currentFace.ownerCell()].push_back
                (
                    currentFace.neighborCell().value()
                );
            }
        }

        // --- Process Neighbor Cell ---
        if 
        (
            currentFace.neighborCell().has_value() 
         && currentFace.neighborCell().value() < allCells.size()
        ) 
        {
            size_t neighborIdx = currentFace.neighborCell().value();
            allCells[neighborIdx].addFaceIndex(faceIdx);
            allCells[neighborIdx].addFaceSign(-1);

            if (currentFace.ownerCell() < allCells.size())
            {
                tempCellNeighbors[neighborIdx].push_back(
                    currentFace.ownerCell());
            }
        }
    }

    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        if (i < allCells.size() && i < tempCellNeighbors.size()) 
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
            
            allCells[i].setNeighborCellIndices(neighbors);
        } 
        else 
        {
            if (i < allCells.size())
            {
                allCells[i].clearNeighborCellIndices();
            }
            std::cerr << "Warning: Cell with index " << i
                      << " could not have its neighbors finalized due to "
                      << "indexing mismatch. allCells.size(): " 
                      << allCells.size()
                      << ", tempCellNeighbors.size(): " 
                      << tempCellNeighbors.size()
                      << ", numCells: " << numCells << std::endl;
        }
    }

    for (size_t i = 0; i < allFaces.size(); ++i) 
    {
        if (allFaces[i].nodeIndices().size() < 3) 
        {
            throw std::runtime_error
            (
                "Error: Face " 
              + std::to_string(i) 
              + " has only " 
              + std::to_string(allFaces[i].nodeIndices().size()) 
              + " nodes, minimum required is 3."
            );
        }
    }
    
    std::cout << "Mesh loaded successfully:" << std::endl;
    std::cout << "  - Nodes: " << allNodes.size() << std::endl;
    std::cout << "  - Faces: " << allFaces.size() << std::endl;
    std::cout << "  - Cells: " << allCells.size() << std::endl;
    std::cout << "  - Boundary patches: " 
              << allBoundaryPatches.size() << std::endl;
}