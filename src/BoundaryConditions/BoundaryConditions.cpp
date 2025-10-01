/******************************************************************************
 * @file BoundaryConditions.cpp
 * @brief Implementation of boundary conditions management system
 *****************************************************************************/

#include "BoundaryConditions.hpp"

#include <stdexcept>
#include <iostream>
#include <algorithm>

void BoundaryConditions::addPatch(const BoundaryPatch& patch)
{
    patches_.push_back(patch);
}

const BoundaryPatch* 
BoundaryConditions::patch(const std::string& name) const
{
    for (const auto& patch : patches_) 
    {
        if (patch.patchName() == name) 
        {
            return &patch;
        }
    }
    
    throw std::runtime_error("Patch " + name + " not found");
}

bool BoundaryConditions::setBC
(
    const std::string& patchName,
    const std::string& fieldName,
    BoundaryData bc_setup
)
{
    patchBoundaryData_[patchName][fieldName] = bc_setup;
    return true;
}

bool BoundaryConditions::setFixedValue
(
    const std::string& patchName,
    const std::string& fieldName,
    Scalar value
)
{
    BoundaryData bc_setup;
    bc_setup.setFixedValue(value);
    return setBC(patchName, fieldName, bc_setup);
}

bool BoundaryConditions::setFixedValue
(
    const std::string& patchName,
    const std::string& fieldName,
    const Vector& value
)
{
    BoundaryData bc_setup;
    bc_setup.setFixedValue(value);
    return setBC(patchName, fieldName, bc_setup);
}

bool BoundaryConditions::setFixedGradient
(
    const std::string& patchName,
    const std::string& fieldName,
    Scalar gradient
)
{
    BoundaryData bc_setup;
    bc_setup.setFixedGradient(gradient);
    return setBC(patchName, fieldName, bc_setup);
}

bool BoundaryConditions::setFixedGradient
(
    const std::string& patchName,
    const std::string& fieldName,
    const Vector& gradient
)
{
    BoundaryData bc_setup;
    bc_setup.setFixedGradient(gradient);
    return setBC(patchName, fieldName, bc_setup);
}

bool BoundaryConditions::setZeroGradient
(
    const std::string& patchName,
    const std::string& fieldName
)
{
    BoundaryData bc_setup;
    bc_setup.setZeroGradient();
    return setBC(patchName, fieldName, bc_setup);
}

bool BoundaryConditions::setNoSlip
(
    const std::string& patchName,
    const std::string& fieldName
)
{
    BoundaryData bc_setup;
    bc_setup.setNoSlip();
    return setBC(patchName, fieldName, bc_setup);
}


const BoundaryData* BoundaryConditions::fieldBC
(
    const std::string& patchName,
    const std::string& fieldName
) const
{
    auto patch_it = patchBoundaryData_.find(patchName);

    if (patch_it != patchBoundaryData_.end()) 
    {
        const auto& field_map = patch_it->second;

        auto field_it = field_map.find(fieldName);

        if (field_it != field_map.end())
        {
            return &(field_it->second);
        }
        
        // If vector component (U_x, U_y, U_z) not found, try parent vector field "U"
        if (fieldName == "U_x" || fieldName == "U_y" || fieldName == "U_z")
        {
            auto vector_field_it = field_map.find("U");
            if (vector_field_it != field_map.end())
            {
                return &(vector_field_it->second);
            }
        }
    }

    return nullptr;
}

Scalar BoundaryConditions::calculateBoundaryFaceValue
(
    const Face& face,
    const ScalarField& phi,
    const std::string& fieldName
) const 
{ 
    // Use shared face-to-patch cache
    ensureFaceToPatchCacheBuilt();
    
    // Find the boundary patch for this face
    auto patch_it = faceToPatchCache_.find(face.id());

    if (patch_it == faceToPatchCache_.end()) 
    {
        throw std::runtime_error
        (
            "Boundary face " 
          + std::to_string(face.id()) 
          + " not found in any boundary patch. Check mesh/BC setup."
        );
    }
    
    const BoundaryPatch* patch = patch_it->second;
    const BoundaryData* bc = fieldBC(patch->patchName(), fieldName);

    if (!bc) {
        // Default to zero-gradient for scalars if not specified
        std::cerr   << "No BC specified for face " << face.id() << " in patch "
                    << patch->patchName() << ". Defaulting to zero-gradient." 
                    << std::endl;

        return phi[face.ownerCell()];
    }
    
    switch (bc->type()) 
    {
        case BCType::FIXED_VALUE:
        {
            // Fixed value: φ_f = φ_b
            return bc->fixedScalarValue();
        }

        case BCType::ZERO_GRADIENT:
        {
            // Zero gradient: φ_f = φ_P
            return phi[face.ownerCell()];
        }
            
        case BCType::FIXED_GRADIENT: 
        {
            // Fixed gradient: φ_f = φ_P + grad * distance
            Scalar d_n = dot(face.d_Pf(), face.normal());
            return phi[face.ownerCell()] + bc->fixedScalarGradient() * d_n;
        }
        
        default:
            throw std::runtime_error
            (
                "Unknown BC type for face " 
              + std::to_string(face.id()) 
              + " in patch " + patch->patchName() 
              + ": " + std::to_string(static_cast<int>(bc->type()))
            );
    }
}

Vector BoundaryConditions::calculateBoundaryFaceVectorValue
(
    const Face& face,
    const VectorField& phi,
    const std::string& fieldName
) const {
        
    // Use shared face-to-patch cache
    ensureFaceToPatchCacheBuilt();
    
    // Find the boundary patch for this face
    auto patch_it = faceToPatchCache_.find(face.id());
    if (patch_it == faceToPatchCache_.end()) 
    {
        throw std::runtime_error
        (
            "Boundary face " 
          + std::to_string(face.id()) 
          + " not found in any boundary patch. Check mesh/BC setup."
        );
    }
    
    const BoundaryPatch* patch = patch_it->second;
    const BoundaryData* bc = fieldBC(patch->patchName(), fieldName);
    if (!bc) {
        // Default to copy owner for vectors if not specified
        std::cerr   << "No BC specified for face " << face.id() << " in patch "
                    << patch->patchName() << ". Defaulting to zero-gradient." 
                    << std::endl;

        return phi[face.ownerCell()];
    }
    
    // Apply boundary condition based on type
    switch (bc->type()) {
        case BCType::FIXED_VALUE:
        {
            // Fixed value: φ_f = φ_b
            return bc->vectorValue();
        }

        case BCType::NO_SLIP:
        {
            // No slip: φ_f = (0.0, 0.0, 0.0)
            return bc->vectorValue(); 
        }

        case BCType::ZERO_GRADIENT:
        {
            // Zero gradient: φ_f = φ_P
            return phi[face.ownerCell()];
        }
            
        case BCType::FIXED_GRADIENT: 
        {
            // Fixed gradient: φ_f = φ_P + grad * distance
            const Vector& d_Pf = face.d_Pf();
            Scalar d_n = dot(d_Pf, face.normal());
            return phi[face.ownerCell()] + bc->vectorGradient() * d_n;
        }

        default:
            throw std::runtime_error
            (
                "Unknown BC type for face " 
              + std::to_string(face.id()) 
              + " in patch " + patch->patchName() 
              + ": " + std::to_string(static_cast<int>(bc->type()))
            );
    }
}

std::string BoundaryConditions::bcTypeToString(BCType bctype) const {
    switch (bctype) {
        case BCType::UNDEFINED: return "UNDEFINED";
        case BCType::FIXED_VALUE: return "FIXED_VALUE";
        case BCType::FIXED_GRADIENT: return "FIXED_GRADIENT";
        case BCType::ZERO_GRADIENT: return "ZERO_GRADIENT";
        case BCType::NO_SLIP: return "NO_SLIP";
        default: 
            throw std::runtime_error
            (
                "Unknown BC type: " 
              + std::to_string(static_cast<int>(bctype))
            );
    }
}

void BoundaryConditions::printSummary() const 
{
    std::cout << "\n--- Boundary Conditions Setup Summary ---" << std::endl;

    if (patches_.empty()) 
    {
        std::cout << "  No mesh patches loaded." << std::endl;
        return;
    }

    std::cout   << "Total Mesh Patches Registered: " << patches_.size() 
                << std::endl;

    for (const auto& meshPatch : patches_) 
    {
        std::cout   << "  ------------------------------------" 
                    << std::endl;

        std::cout   << "  Mesh Patch Name         : " 
                    << meshPatch.patchName() << std::endl;

        std::cout   << "  Fluent Type             : " 
                    << meshPatch.fluentType() << std::endl;

        std::cout   << "  Zone ID                 : " 
                    << meshPatch.zoneID() << std::endl;

        std::cout   << "  Number of Faces         : " 
                    << meshPatch.numberOfBoundaryFaces() << std::endl;

        auto patch_bc_it = patchBoundaryData_.find(meshPatch.patchName());

        if 
        (
            patch_bc_it != patchBoundaryData_.end() 
         && !patch_bc_it->second.empty()
        ) 
        {
            std::cout << "  Configured Physical BCs :" << std::endl;

            for (const auto& field_bc_pair : patch_bc_it->second) 
            {
                const std::string& fieldName = field_bc_pair.first;
                const BoundaryData& fbc = field_bc_pair.second;
                
                std::cout   << "      Field '" << fieldName << "': Type: "
                            << bcTypeToString(fbc.type());
                
                if 
                (
                    fbc.type() == BCType::FIXED_VALUE 
                 || fbc.type() == BCType::NO_SLIP
                )
                {
                    std::cout << ", Value: ";

                    if (fbc.valueType() == BCValueType::SCALAR)
                    {
                        std::cout << fbc.scalarValue();
                    }
                    else if
                    (fbc.valueType() == BCValueType::VECTOR) {
                        std::cout << fbc.vectorValue();
                    }
                    else
                    {
                        throw std::runtime_error
                        (
                            "Unknown BC value type: " 
                          + std::to_string(static_cast<int>(fbc.valueType()))
                        );
                    }
                } 
                else if (fbc.type() == BCType::FIXED_GRADIENT)
                {
                    std::cout << ", Gradient: ";

                    if (fbc.gradientType() == BCValueType::SCALAR)
                    {
                        std::cout << fbc.scalarGradient();
                    } 
                    else if (fbc.gradientType() == BCValueType::VECTOR)
                    {
                        std::cout << fbc.vectorGradient();
                    } 
                    else
                    {
                        throw std::runtime_error
                        (
                            "Unknown BC gradient type: "
                          + std::to_string(static_cast<int>(fbc.gradientType()))
                        );
                    }
                }
                else if (fbc.type() == BCType::ZERO_GRADIENT)
                {
                    std::cout << " (implies zero gradient)";
                }
                std::cout << std::endl;
            }
        }
    }

    std::cout << "  ------------------------------------" << std::endl;
}

void BoundaryConditions::ensureFaceToPatchCacheBuilt() const
{
    if (!cacheBuilt_)
    {
        faceToPatchCache_.clear();
        
        for (const auto& patch : patches_)
        {
            for 
            (
                size_t f = patch.firstFaceIndex(); f <= patch.lastFaceIndex(); ++f
            )
            {
                faceToPatchCache_[f] = &patch;
            }
        }

        cacheBuilt_ = true;
    }
}
