#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>    // For std::numeric_limits
#include <algorithm> // For std::sort, std::unique
#include <stdexcept>
#include <charconv>

#include "MeshReader.h"

// Section codes in Fluent file 
const std::string MSH_COMMENT = "(0";
const std::string MSH_DIMENSION = "(2";
const std::string MSH_NODES = "(10";
const std::string MSH_CELLS = "(12";
const std::string MSH_FACES = "(13";
const std::string MSH_BOUNDARIES = "(45";

// Helper to convert hex string to size_t
size_t hexToDec(const std::string &hexStr) {
    size_t decVal;
    auto result = std::from_chars(hexStr.data(), hexStr.data() + hexStr.size(), decVal, 16);
    
    // Error checking: ensure the entire string was parsed
    if (result.ec != std::errc() || result.ptr != hexStr.data() + hexStr.size()) {
        throw std::runtime_error("Error: Failed to convert hex string '" + hexStr + "' to decimal.");
    }
    return decVal;
}

// Helper to convert decimal string to size_t
size_t strToDec(const std::string& decStr) {
    if (decStr.empty()) {
        throw std::runtime_error("Error: Attempted to convert empty decimal string to size_t.");
    }
    std::stringstream ss;
    ss << decStr; 
    size_t decVal;
    ss >> decVal;
    if (ss.fail() || !ss.eof()) {

        throw std::runtime_error("Error: Failed to convert decimal string '" + decStr + "' to size_t. Invalid format or contains non-numeric characters.");
    }
    return decVal;
}

void readMshFile(const std::string &filePath,
                 std::vector<Vector> &allNodes,
                 std::vector<Face> &allFaces,
                 std::vector<Cell> &allCells,
                 std::vector<BoundaryPatch> &allBoundaryPatches)
{

  // Ensure vectors are empty before reading
  allNodes.clear();
  allFaces.clear();
  allCells.clear();
  allBoundaryPatches.clear();

  std::string line;
  std::string token;
  char ch;

  size_t num_nodes = 0;
  size_t num_faces = 0;
  size_t num_cells = 0;

  bool is2D = false;

  std::ifstream ifs(filePath);
  if (!ifs.is_open()) {
    throw std::runtime_error("Error: Could not open mesh file: " + filePath);
  }

  while (ifs >> token) {
    // ----- Comment Section ----- //
    if (token == MSH_COMMENT) { 
      int paren_level = 1;

      // Switch from token to character 
      while (ifs.get(ch)) {
        if (ch == '(')
          paren_level++;
        else if (ch == ')')
          paren_level--;
        if (paren_level == 0)
          break;
      }
    }

    // ----- Dimension Section ----- //
    if (token == MSH_DIMENSION) { 
      std::string dimension;
      ifs >> dimension;
      dimension.pop_back();
      if (dimension == "2") {
        is2D = true;
        throw std::runtime_error("This code doesn't handle 2D geometries!");
      } else if (dimension == "3")
        is2D = false;
    }

    // ----- Points Section ----- //
    if (token == MSH_NODES) { 
      ifs >> token;

      // Declaration at the beginning of the file 
      if (token == "(0") {
        std::string first_idx_str, last_idx_str, zero_str_or_type;
        ifs >> first_idx_str;
        ifs >> last_idx_str;
        ifs >> zero_str_or_type;

        size_t last_idx = hexToDec(last_idx_str);

        allNodes.resize(last_idx);
        num_nodes = allNodes.size();
      }

      // Points data
      else {
        std::string zone_id_str, start_id_str, end_id_str, type_str, dimension_str;
        
        token.erase(0, 1);
        zone_id_str = token;

        ifs >> start_id_str;
        ifs >> end_id_str;
        ifs >> type_str;
        ifs >> dimension_str;

        dimension_str.pop_back();
        dimension_str.pop_back();

        size_t start_id = hexToDec(start_id_str);
        size_t end_id = hexToDec(end_id_str);
        size_t dimension = hexToDec(dimension_str);

        size_t node_global_idx = start_id - 1; // 0-based

        for (size_t i = start_id; i <= end_id; ++i) {
          if (node_global_idx < end_id) {
            ifs >> allNodes[node_global_idx].x;
            ifs >> allNodes[node_global_idx].y;
            if (dimension == 3) {
              ifs >> allNodes[node_global_idx].z;
            }
            node_global_idx++;
          } else {
            throw std::runtime_error("Error: Node index " + std::to_string(node_global_idx) + 
                                   " exceeds allocated node vector size " + std::to_string(allNodes.size()));            
          }
        }
      }
    }

    // ----- Cells Section ----- //
    // Cells are not defined explicitly in the file. However, defined by data in the Faces section  
    if (token == MSH_CELLS) {
      ifs >> token;
      
      // Declaration at the beginning of the file 
      if (token == "(0") {
        std::string first_idx_str, last_idx_str, zero_str_or_type;
        ifs >> first_idx_str;
        ifs >> last_idx_str;
        ifs >> zero_str_or_type;

        size_t last_idx = hexToDec(last_idx_str);

        allCells.resize(last_idx);
        num_cells = allCells.size();
      }
    }

    // ----- Faces Section ----- //
    if (token == MSH_FACES) {
      ifs >> token;

      // Declaration at the beginning of the file 
      if (token == "(0") {
        std::string first_idx_str, last_idx_str, zero_str_or_type;
        ifs >> first_idx_str;
        ifs >> last_idx_str;
        ifs >> zero_str_or_type;

        size_t last_idx = hexToDec(last_idx_str);

        allFaces.resize(last_idx);
        num_faces = allFaces.size();
      }

      // Faces data
      else {
        std::string zone_id_str, start_id_str, end_id_str, type_str, elementType_str;
        token.erase(0, 1);
        zone_id_str = token;

        ifs >> start_id_str;
        ifs >> end_id_str;
        ifs >> type_str;
        ifs >> elementType_str;

        elementType_str.pop_back();
        elementType_str.pop_back();

        size_t zone_id = hexToDec(zone_id_str);
        size_t start_id = hexToDec(start_id_str);
        size_t end_id = hexToDec(end_id_str);
        
        allBoundaryPatches.emplace_back(zone_id, start_id-1, end_id-1);

        // Switch from token to line 
        std::string line;
        std::getline(ifs, line); 

        for (size_t face_i = start_id; face_i <= end_id; ++face_i) {
          if ((face_i - 1) >= allFaces.size()) {
            throw std::runtime_error("Error: Face index " + std::to_string((face_i - 1)) +
                          " is out of bounds for allFaces vector of size " + std::to_string(allFaces.size()) + ".");
          }

          Face &currentFace = allFaces[(face_i - 1)];
          currentFace.id = face_i - 1;
          currentFace.nodeIndices.clear();

          if (!std::getline(ifs, line)) {
            throw std::runtime_error("Error: Could not read face data line for face ID (1-based): " + std::to_string(face_i) + 
                                   ". End of file reached unexpectedly.");
          }
          
          if (line.empty()) {
            throw std::runtime_error("Error: Empty line encountered while reading face data for face ID (1-based): " + std::to_string(face_i));
          }
          
          std::stringstream line_stream(line);
          std::string item_hex;
          std::vector<std::string> hex_items;

          // Put hex items in a vector 
          while (line_stream >> item_hex) {
            hex_items.push_back(item_hex);
          }

          // The last two hex ids are for the owner cell and the neighbor cell.
          std::string neighbour_hex = hex_items.back();
          hex_items.pop_back();
          std::string owner_hex = hex_items.back();
          hex_items.pop_back();

          // The remaining items are node indices
          for (const std::string &node_idx_hex : hex_items) {
            currentFace.nodeIndices.push_back(hexToDec(node_idx_hex) - 1);
          }

          currentFace.ownerCell = hexToDec(owner_hex) - 1;

          if (neighbour_hex != "0") {
            currentFace.neighbourCell = hexToDec(neighbour_hex) - 1;
          }
          else {
            currentFace.neighbourCell = std::nullopt; // Boundary
          }
        } 
      }
    }

    // ----- Boundaries Section ----- //
    if (token == MSH_BOUNDARIES) {
      ifs >> token;

      // Remove leading '('
      token.erase(0, 1); 
      size_t zoneID = strToDec(token);

      std::string typeString;
      ifs >> typeString;

      std::string nameString;
      ifs >> nameString;
      nameString.pop_back();
      nameString.pop_back();
      nameString.pop_back();
      nameString.pop_back();

      for (size_t i = 0; i < allBoundaryPatches.size(); i++) {
        if (zoneID == allBoundaryPatches[i].zoneID) {
          allBoundaryPatches[i].patchName = nameString;
          allBoundaryPatches[i].fluentType = typeString;
          BoundaryConditionType mappedType = mapFluentBCToEnum(typeString);
          allBoundaryPatches[i].type = mappedType;
        }
      }
    }
  }

  ifs.close();

  // ----- Populate Cell Data (Face Indices and Neighbour Cell Indices) ----- //

  // Initialize cell's ID and pre-allocate vectors for better performance
  for (size_t i = 0; i < num_cells; ++i) {
    if (i < allCells.size())
    {
      allCells[i].id = i;
      allCells[i].faceIndices.clear();
      allCells[i].neighbourCellIndices.clear();
    }
    else {
      throw std::runtime_error("Error: Cell vector out of bounds during initialization. num_cells: " + std::to_string(num_cells) + ", index: " + std::to_string(i));
    }
  }

  // Temporary storage for accumulating unique neighbor cell indices for each cell.
  std::vector<std::vector<size_t>> tempCellNeighbors(num_cells);

  // Iterate through all faces to establish cell-to-face and cell-to-cell relationships.
  for (size_t face_idx = 0; face_idx < allFaces.size(); ++face_idx) {
    const Face &currentFace = allFaces[face_idx];

    // --- Process Owner Cell ---
    if (currentFace.ownerCell != std::numeric_limits<size_t>::max() && currentFace.ownerCell < num_cells) {
      // Add this face's index to its owner cell's list of faces.
      allCells[currentFace.ownerCell].faceIndices.push_back(face_idx);
      allCells[currentFace.ownerCell].faceSigns.push_back(1);

      // If this face has a valid neighbourCell, then that neighbourCell
      if (currentFace.neighbourCell.has_value() && currentFace.neighbourCell.value() < num_cells) {
        tempCellNeighbors[currentFace.ownerCell].push_back(currentFace.neighbourCell.value());
      }
    }

    // --- Process Neighbour Cell ---
    if (currentFace.neighbourCell.has_value() && currentFace.neighbourCell.value() < num_cells) {
      size_t neighborIdx = currentFace.neighbourCell.value();
      allCells[neighborIdx].faceIndices.push_back(face_idx);
      allCells[neighborIdx].faceSigns.push_back(-1);

      // The ownerCell of this face is a neighboring cell to the current 'neighbourCell'.
      // This connection is only made if the ownerCell itself is a valid internal cell.
      if (currentFace.ownerCell < num_cells)
      {
        tempCellNeighbors[neighborIdx].push_back(currentFace.ownerCell);
      }
    }
  }

  // Sort and remove duplicate neighbor indices for each cell.
  for (size_t i = 0; i < num_cells; ++i) {
    if (i < allCells.size() && i < tempCellNeighbors.size()) {
      std::vector<size_t> &neighbors = tempCellNeighbors[i];
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
      allCells[i].neighbourCellIndices = neighbors; // Assign the processed list
    } else {
      if (i < allCells.size())
      {
        allCells[i].neighbourCellIndices.clear();
      }
      std::cerr << "Warning: Cell with index " << i
                << " could not have its neighbors finalized due to indexing mismatch. "
                << "allCells.size(): " << allCells.size()
                << ", tempCellNeighbors.size(): " << tempCellNeighbors.size()
                << ", num_cells: " << num_cells << std::endl;
    }
  }

  // ----- Final Mesh Validation ----- //
  
  // Validate that all cells have at least 4 faces (minimum for 3D)
  for (size_t i = 0; i < allCells.size(); ++i) {
    if (allCells[i].faceIndices.size() < 4) {
      std::cerr << "Warning: Cell " << i << " has only " << allCells[i].faceIndices.size() 
                << " faces, which is below the minimum of 4 for 3D cells." << std::endl;
    }
  }
  
  // Validate that all faces have at least 3 nodes
  for (size_t i = 0; i < allFaces.size(); ++i) {
    if (allFaces[i].nodeIndices.size() < 3) {
      throw std::runtime_error("Error: Face " + std::to_string(i) + " has only " + 
                             std::to_string(allFaces[i].nodeIndices.size()) + " nodes, minimum required is 3.");
    }
  }
  
  // Print summary statistics
  std::cout << "Mesh loaded successfully:" << std::endl;
  std::cout << "  - Nodes: " << allNodes.size() << std::endl;
  std::cout << "  - Faces: " << allFaces.size() << std::endl;
  std::cout << "  - Cells: " << allCells.size() << std::endl;
  std::cout << "  - Boundary patches: " << allBoundaryPatches.size() << std::endl;
}