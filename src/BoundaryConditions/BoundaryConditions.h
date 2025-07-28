#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <vector>
#include <string>
#include <map>
#include <stdexcept>
#include <iostream>

#include "BoundaryPatch.h"
#include "BoundaryData.h"
#include "Vector.h"
#include "Scalar.h"

class BoundaryConditions {
public:
    std::vector<BoundaryPatch> patches;

    // The main storage of boundary conditions setup 
    std::map<std::string, std::map<std::string, BoundaryData>> patchBoundaryData;

    // ----- Constructors ------//
    BoundaryConditions() = default;

    // Add patch from the mesh reader
    void addPatch(const BoundaryPatch& patch) {
        patches.push_back(patch);
    }

    // Get a patch by its name
    const BoundaryPatch* getPatch(const std::string& name) const {
        for (const auto& patch : patches) {
            if (patch.patchName == name) {
                return &patch;
            }
        }
        return nullptr; // Not found
    }


    // ----- Methods to set boundary conditions ----- //

    bool setBC(const std::string& patchName, const std::string& fieldName, BoundaryData bc_config) {
        if (getPatch(patchName) == nullptr) {
            return false;
        }
        patchBoundaryData[patchName][fieldName] = bc_config;
        return true;
    }

    bool setFixedValue(const std::string& patchName, const std::string& fieldName, Scalar value) {
        BoundaryData bc_config;
        bc_config.setFixedValue(value);
        return setBC(patchName, fieldName, bc_config);
    }

    bool setFixedValue(const std::string& patchName, const std::string& fieldName, const Vector& value) {
        BoundaryData bc_config;
        bc_config.setFixedValue(value);
        return setBC(patchName, fieldName, bc_config);
    }

    bool setZeroGradient(const std::string& patchName, const std::string& fieldName) {
        BoundaryData bc_config;
        bc_config.setZeroGradient();
        return setBC(patchName, fieldName, bc_config);
    }
    
    bool setNoSlip(const std::string& patchName, const std::string& fieldName) {
        BoundaryData bc_config;
        bc_config.setNoSlip();
        return setBC(patchName, fieldName, bc_config);
    }

    bool setFixedGradient(const std::string& patchName, const std::string& fieldName, Scalar gradient) {
        BoundaryData bc_config;
        bc_config.setFixedGradient(gradient);
        return setBC(patchName, fieldName, bc_config);
    }

    bool setFixedGradient(const std::string& patchName, const std::string& fieldName, const Vector& gradient) {
        BoundaryData bc_config;
        bc_config.setFixedGradient(gradient);
        return setBC(patchName, fieldName, bc_config);
    }

    // ----- Methods to get boundary condition configuration ----- //

    const BoundaryData* getFieldBC(const std::string& patchName, const std::string& fieldName) const {
        auto patch_it = patchBoundaryData.find(patchName);
        if (patch_it != patchBoundaryData.end()) {
            const auto& field_map = patch_it->second;
            auto field_it = field_map.find(fieldName);
            if (field_it != field_map.end()) {
                return &(field_it->second);
            }
        }
        return nullptr; 
    }

    // Helper Method to convert BCType enum to string for printing
    std::string bcTypeToString(BCType bctype) const {
        switch (bctype) {
            case BCType::UNDEFINED: return "UNDEFINED";
            case BCType::FIXED_VALUE: return "FIXED_VALUE";
            case BCType::FIXED_GRADIENT: return "FIXED_GRADIENT";
            case BCType::ZERO_GRADIENT: return "ZERO_GRADIENT";
            case BCType::NO_SLIP: return "NO_SLIP";
            default: return "UNKNOWN_BCTYPE";
        }
    }

    // Method to print boundary conditions data 
    void printSummary(bool detailed = true) const {
        std::cout << "\n--- Boundary Conditions Setup Summary ---" << std::endl;
        if (patches.empty()) {
            std::cout << "  No mesh patches loaded." << std::endl;
            return;
        }
        std::cout << "Total Mesh Patches Registered: " << patches.size() << std::endl;
        for (const auto& meshPatch : patches) {
            std::cout << "  ------------------------------------" << std::endl;
            std::cout << "  Mesh Patch Name    : " << meshPatch.patchName << std::endl;
            std::cout << "  Fluent Type        : " << meshPatch.fluentType << std::endl;
            std::cout << "  Zone ID            : " << meshPatch.zoneID << std::endl;
            
            if (detailed) {
                if (meshPatch.firstFaceIndex != std::numeric_limits<size_t>::max()) {
                    std::cout << "  File Block Range   : [" << meshPatch.firstFaceIndex
                              << " - " << meshPatch.lastFaceIndex << "]" << std::endl;
                } else {
                     std::cout << "  File Block Range   : Not set" << std::endl;
                }
                std::cout << "  Total Faces in Block: " << meshPatch.getNumberOfBoundaryFaces() << std::endl;

                size_t numFacesInPatch = meshPatch.getNumberOfBoundaryFaces();
                if (numFacesInPatch > 0) {
                    std::cout << "    First few boundary face IDs: [";
                    size_t numToPrint = std::min(numFacesInPatch, (size_t)5);
                    for (size_t i = 0; i < numToPrint; ++i) {
                        // Print the face ID by adding the offset to the starting index
                        std::cout << (meshPatch.firstFaceIndex + i);
                        if (i < numToPrint - 1) {
                            std::cout << ", ";
                        }
                    }
                    if (numFacesInPatch > 5) {
                        std::cout << ", ...";
                    }
                    std::cout << "]" << std::endl;
                }
            }

            // Print the configured physical BCs for this patch
            auto patch_bc_it = patchBoundaryData.find(meshPatch.patchName);
            if (patch_bc_it != patchBoundaryData.end() && !patch_bc_it->second.empty()) {
                std::cout << "    Configured Physical BCs:" << std::endl;
                for (const auto& field_bc_pair : patch_bc_it->second) {
                    const std::string& fieldName = field_bc_pair.first;
                    const BoundaryData& fbc = field_bc_pair.second;
                    
                    std::cout << "      Field '" << fieldName << "': Type: " << bcTypeToString(fbc.type);
                    
                    if (fbc.type == BCType::FIXED_VALUE || fbc.type == BCType::NO_SLIP) {
                        std::cout << ", Value: ";
                        if (fbc.valueType == BCValueType::SCALAR) {
                            std::cout << fbc.scalarValue;
                        } else if (fbc.valueType == BCValueType::VECTOR) {
                            std::cout << fbc.vectorValue;
                        } else {
                            std::cout << "N/A (value type not set)";
                        }
                    } else if (fbc.type == BCType::FIXED_GRADIENT) {
                        std::cout << ", Gradient: ";
                        if (fbc.gradientType == BCValueType::SCALAR) {
                            std::cout << fbc.scalarGradient;
                        } else if (fbc.gradientType == BCValueType::VECTOR) {
                            std::cout << fbc.vectorGradient;
                        } else {
                            std::cout << "N/A (gradient type not set)";
                        }
                    } else if (fbc.type == BCType::ZERO_GRADIENT) {
                        // No specific value to print for ZERO_GRADIENT itself beyond its type
                        std::cout << " (implies zero gradient)";
                    }
                    std::cout << std::endl;
                }
            } else {
                 std::cout << "    No specific physical BCs configured for this patch name." << std::endl;
            }
        }
        std::cout << "  ------------------------------------" << std::endl;
    }
};

#endif
