#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <vector>
#include <string>
#include <map>
#include <stdexcept>
#include <iostream>

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"
#include "BoundaryData.h"
#include "CellData.h"

class BoundaryConditions {
public:
    // ----- Constructor ----- //
    BoundaryConditions() = default;

    // ----- Methods ----- //

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
        // if not found
        throw std::runtime_error("Patch " + name + " not found");
    }

    // ----- Boundary Condition Setters ----- //
    // Generic setter
    bool setBC(const std::string& patchName, const std::string& fieldName, BoundaryData bc_config) {
        patchBoundaryData[patchName][fieldName] = bc_config;
        return true;
    }

    // Fixed value setters
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

    // Gradient setters
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

    // Special boundary condition setters
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

    bool setSymmetry(const std::string& patchName, const std::string& fieldName) {
        BoundaryData bc_config;
        bc_config.setSymmetry();
        return setBC(patchName, fieldName, bc_config);
    }

    // ----- Boundary Condition Getters ----- //
    // Get boundary condition configuration
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
        ;
    }

    // ----- Boundary Face Value Calculations ----- //
    // Calculate boundary face value based on boundary conditions
    Scalar calculateBoundaryFaceValue(
        const Face& face,
        const ScalarField& phi,
        const std::string& fieldName) const {
            
        // Use shared face-to-patch cache
        ensureFaceToPatchCacheBuilt();
        
        // Find the boundary patch for this face
        auto patch_it = faceToPatchCache.find(face.id);
        if (patch_it == faceToPatchCache.end()) {
            throw std::runtime_error("Boundary face " + std::to_string(face.id) + 
                                   " not found in any boundary patch. Check mesh/BC setup.");
        }
        
        const BoundaryPatch* patch = patch_it->second;
        const BoundaryData* bc = getFieldBC(patch->patchName, fieldName);
        if (!bc) {
            // Default to zero-gradient for scalars if not specified
            return phi[face.ownerCell];
        }
        
        // Apply boundary condition based on type
        switch (bc->type) {
            case BCType::FIXED_VALUE:
                return bc->getFixedScalarValue();

            case BCType::ZERO_GRADIENT:
                // Zero gradient: φ_f = φ_P
                return phi[face.ownerCell];
                
            case BCType::FIXED_GRADIENT: {
                // Fixed gradient: φ_f = φ_P + grad * distance
                Scalar d_n = dot(face.d_Pf, face.normal);
                return phi[face.ownerCell] + bc->getFixedScalarGradient() * d_n;
            }
            case BCType::SYMMETRY:
                // Symmetry: zero normal gradient for scalars
                return phi[face.ownerCell];
            
            default:
                throw std::runtime_error("Unknown BC type for face " + std::to_string(face.id) + 
                                       " in patch " + patch->patchName + ": " + std::to_string(static_cast<int>(bc->type)));
        }
    }

    // Calculate boundary face vector value based on boundary conditions
    Vector calculateBoundaryFaceVectorValue(
        const Face& face,
        const VectorField& phi,
        const std::string& fieldName) const {
            
        // Use shared face-to-patch cache
        ensureFaceToPatchCacheBuilt();
        
        // Find the boundary patch for this face
        auto patch_it = faceToPatchCache.find(face.id);
        if (patch_it == faceToPatchCache.end()) {
            throw std::runtime_error("Boundary face " + std::to_string(face.id) + 
                                   " not found in any boundary patch. Check mesh/BC setup.");
        }
        
        const BoundaryPatch* patch = patch_it->second;
        const BoundaryData* bc = getFieldBC(patch->patchName, fieldName);
        if (!bc) {
            // Default to copy owner for vectors if not specified
            return phi[face.ownerCell];
        }
        
        // Apply boundary condition based on type
        switch (bc->type) {
            case BCType::FIXED_VALUE:
                return bc->vectorValue;
                
            case BCType::NO_SLIP:
                return bc->vectorValue; 
                
            case BCType::ZERO_GRADIENT:
                // Zero gradient: φ_f = φ_P
                return phi[face.ownerCell];
            case BCType::SYMMETRY: {
                // Zero normal component at symmetry plane: U_f = U_P - (U_P·n) n
                Vector U_P = phi[face.ownerCell];
                Vector n = face.normal;
                return U_P - dot(U_P, n) * n;
            }
                
            case BCType::FIXED_GRADIENT: {
                // Fixed gradient: φ_f = φ_P + grad * distance
                const Vector& d_Pf = face.d_Pf;
                Scalar d_n = dot(d_Pf, face.normal);
                return phi[face.ownerCell] + bc->vectorGradient * d_n;
            }
            
            default:
                throw std::runtime_error("Unknown BC type for face " + std::to_string(face.id) + 
                                       " in patch " + patch->patchName + ": " + std::to_string(static_cast<int>(bc->type)));
        }
    }

    // ----- Utility & Debugging Methods ----- //
    // Convert BCType enum to string for printing
    std::string bcTypeToString(BCType bctype) const {
        switch (bctype) {
            case BCType::UNDEFINED: return "UNDEFINED";
            case BCType::FIXED_VALUE: return "FIXED_VALUE";
            case BCType::FIXED_GRADIENT: return "FIXED_GRADIENT";
            case BCType::ZERO_GRADIENT: return "ZERO_GRADIENT";
            case BCType::NO_SLIP: return "NO_SLIP";
            case BCType::SYMMETRY: return "SYMMETRY";
            default: 
                throw std::runtime_error("Unknown BC type: " + std::to_string(static_cast<int>(bctype)));
        }
    }

    // Print boundary conditions data 
    void printSummary() const {
        std::cout << "\n--- Boundary Conditions Setup Summary ---" << std::endl;
        if (patches.empty()) {
            std::cout << "  No mesh patches loaded." << std::endl;
            return;
        }
        std::cout << "Total Mesh Patches Registered: " << patches.size() << std::endl;
        for (const auto& meshPatch : patches) {
            std::cout << "  ------------------------------------" << std::endl;
            std::cout << "  Mesh Patch Name         : " << meshPatch.patchName << std::endl;
            std::cout << "  Fluent Type             : " << meshPatch.fluentType << std::endl;
            std::cout << "  Zone ID                 : " << meshPatch.zoneID << std::endl;
            std::cout << "  Number of Faces         : " << meshPatch.getNumberOfBoundaryFaces() << std::endl;

            // Print the configured physical BCs for this patch
            auto patch_bc_it = patchBoundaryData.find(meshPatch.patchName);
            if (patch_bc_it != patchBoundaryData.end() && !patch_bc_it->second.empty()) {
                std::cout << "  Configured Physical BCs :" << std::endl;
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
                            throw std::runtime_error("Unknown BC value type: " + std::to_string(static_cast<int>(fbc.valueType)));
                        }
                    } else if (fbc.type == BCType::FIXED_GRADIENT) {
                        std::cout << ", Gradient: ";
                        if (fbc.gradientType == BCValueType::SCALAR) {
                            std::cout << fbc.scalarGradient;
                        } else if (fbc.gradientType == BCValueType::VECTOR) {
                            std::cout << fbc.vectorGradient;
                        } else {
                            throw std::runtime_error("Unknown BC gradient type: " + std::to_string(static_cast<int>(fbc.gradientType)));
                        }
                    } else if (fbc.type == BCType::ZERO_GRADIENT) {
                        std::cout << " (implies zero gradient)";
                    }
                    std::cout << std::endl;
                }
            }
        }
        std::cout << "  ------------------------------------" << std::endl;
    }

private:
    // ----- Private Member Variables ----- //
    // Shared cache for face-to-patch mapping
    mutable std::map<size_t, const BoundaryPatch*> faceToPatchCache;
    mutable bool cacheBuilt = false;
    
    // ----- Private Helper Methods ----- //
    // Helper method to ensure cache is built
    void ensureFaceToPatchCacheBuilt() const {
        if (!cacheBuilt) {
            faceToPatchCache.clear();
            for (const auto& patch : patches) {
                for (size_t f = patch.firstFaceIndex; f <= patch.lastFaceIndex; ++f) {
                    faceToPatchCache[f] = &patch;
                }
            }
            cacheBuilt = true;
        }
    }

public:
    // ----- Public Member Variables ----- //
    // The main storage of boundary conditions setup
    /* Example:
    bcManager
    ├── "inlet"
    │   ├── "U"     --> BoundaryData(FIXED_VALUE, VECTOR, (10,0,0))
    │   ├── "p"     --> BoundaryData(ZERO_GRADIENT)
    │   ├── "k"     --> BoundaryData(FIXED_VALUE, SCALAR, 0.1)
    │   └── "omega" --> BoundaryData(FIXED_VALUE, SCALAR, 1.0)
    │
    ├── "outlet"
    │   ├── "U"     --> BoundaryData(ZERO_GRADIENT)
    │   ├── "p"     --> BoundaryData(FIXED_VALUE, SCALAR, 101325.0)
    │   ├── "k"     --> BoundaryData(ZERO_GRADIENT)
    │   └── "omega" --> BoundaryData(ZERO_GRADIENT)
    */
    std::map<std::string, std::map<std::string, BoundaryData>> patchBoundaryData;
    std::vector<BoundaryPatch> patches;
};

#endif