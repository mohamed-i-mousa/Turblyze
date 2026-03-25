/******************************************************************************
 * @file BoundaryConditions.cpp
 * @brief Implementation of boundary conditions management system
 *****************************************************************************/

#include "BoundaryConditions.hpp"

#include <stdexcept>
#include <iostream>
#include <utility>
#include <set>


// ****************************** Setter Methods ******************************

void BoundaryConditions::addPatch(BoundaryPatch patch)
{
    patches_.push_back(std::move(patch));
}

bool BoundaryConditions::setBC
(
    const std::string& patchName,
    const std::string& fieldName,
    BoundaryData bcData
)
{
    patchBoundaryData_[patchName][fieldName] = std::move(bcData);
    return true;
}

bool BoundaryConditions::setFixedValue
(
    const std::string& patchName,
    const std::string& fieldName,
    Scalar value
)
{
    BoundaryData bcData;
    bcData.setFixedValue(value);
    return setBC(patchName, fieldName, std::move(bcData));
}

bool BoundaryConditions::setFixedValue
(
    const std::string& patchName,
    const std::string& fieldName,
    const Vector& value
)
{
    BoundaryData bcData;
    bcData.setFixedValue(value);
    return setBC(patchName, fieldName, std::move(bcData));
}

bool BoundaryConditions::setFixedGradient
(
    const std::string& patchName,
    const std::string& fieldName,
    Scalar gradient
)
{
    BoundaryData bcData;
    bcData.setFixedGradient(gradient);
    return setBC(patchName, fieldName, std::move(bcData));
}

bool BoundaryConditions::setFixedGradient
(
    const std::string& patchName,
    const std::string& fieldName,
    const Vector& gradient
)
{
    BoundaryData bcData;
    bcData.setFixedGradient(gradient);
    return setBC(patchName, fieldName, std::move(bcData));
}

bool BoundaryConditions::setZeroGradient
(
    const std::string& patchName,
    const std::string& fieldName
)
{
    BoundaryData bcData;
    bcData.setZeroGradient();
    return setBC(patchName, fieldName, std::move(bcData));
}

bool BoundaryConditions::setNoSlip
(
    const std::string& patchName,
    const std::string& fieldName
)
{
    BoundaryData bcData;
    bcData.setNoSlip();
    return setBC(patchName, fieldName, std::move(bcData));
}

bool BoundaryConditions::setKWallFunction
(
    const std::string& patchName,
    const std::string& fieldName
)
{
    BoundaryData bcData;
    bcData.setKWallFunction();
    return setBC(patchName, fieldName, std::move(bcData));
}

bool BoundaryConditions::setOmegaWallFunction
(
    const std::string& patchName,
    const std::string& fieldName
)
{
    BoundaryData bcData;
    bcData.setOmegaWallFunction();
    return setBC(patchName, fieldName, std::move(bcData));
}

bool BoundaryConditions::setNutWallFunction
(
    const std::string& patchName,
    const std::string& fieldName
)
{
    BoundaryData bcData;
    bcData.setNutWallFunction();
    return setBC(patchName, fieldName, std::move(bcData));
}


// ****************************** Accessor Methods ******************************

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

    throw
        std::runtime_error
        (
            "Patch " + name + " not found"
        );
}

const BoundaryData* BoundaryConditions::fieldBC
(
    const std::string& patchName,
    const std::string& fieldName
) const
{
    auto patchIterator = patchBoundaryData_.find(patchName);

    if (patchIterator != patchBoundaryData_.end())
    {
        const auto& fieldMap = patchIterator->second;

        auto fieldIterator = fieldMap.find(fieldName);

        if (fieldIterator != fieldMap.end())
        {
            return &(fieldIterator->second);
        }
    }

    return nullptr;
}

Scalar BoundaryConditions::calculateBoundaryFaceValue
(
    const std::string& fieldName,
    const ScalarField& phi,
    const Face& face,
    std::optional<int> componentIdx
) const
{
    const BoundaryPatch* patch = face.patch();
    const BoundaryData* bc = fieldBC(patch->patchName(), fieldName);

    if (!bc)
    {
        std::cerr
            << "No BC specified for face "
            << face.idx() << " in patch "
            << patch->patchName()
            << ". Defaulting to zero-gradient."
            << std::endl;

        return phi[face.ownerCell()];
    }

    using enum BCType;
    switch (bc->type())
    {
        case NO_SLIP:
        case FIXED_VALUE:
        {
            // Extract scalar from vector BC using component index
            if (componentIdx && bc->valueType() == BCValueType::VECTOR)
            {
                const Vector& v = bc->fixedVectorValue();
                switch (*componentIdx)
                {
                    case 0: return v.x();
                    case 1: return v.y();
                    case 2: return v.z();
                }
            }

            // Fixed value: φf = φb
            return bc->fixedScalarValue();
        }

        case K_WALL_FUNCTION:
        case OMEGA_WALL_FUNCTION:
        case NUT_WALL_FUNCTION:
        case ZERO_GRADIENT:
        {
            // Zero gradient: φf = φP
            return phi[face.ownerCell()];
        }

        case FIXED_GRADIENT:
        {
            // Fixed gradient: φf = φP + grad * distance
            Scalar dn = dot(face.dPf(), face.normal());
            return
                phi[face.ownerCell()]
              + bc->fixedScalarGradient() * dn;
        }

        default:
            throw
                std::runtime_error
                (
                    "Unknown BC type for face "
                  + std::to_string(face.idx())
                  + " in patch " + patch->patchName()
                  + ": "
                  + std::to_string(static_cast<int>(bc->type()))
                );
    }
}

Vector BoundaryConditions::calculateBoundaryVectorFaceValue
(
    const std::string& fieldName,
    const VectorField& phi,
    const Face& face
) const
{
    const BoundaryPatch* patch = face.patch();
    const BoundaryData* bc = fieldBC(patch->patchName(), fieldName);

    if (!bc)
    {
        std::cerr
            << "No BC specified for face "
            << face.idx() << " in patch "
            << patch->patchName()
            << ". Defaulting to zero-gradient."
            << std::endl;

        return phi[face.ownerCell()];
    }

    using enum BCType;
    switch (bc->type())
    {
        case NO_SLIP:
        {
            return Vector(0.0, 0.0, 0.0);
        }

        case FIXED_VALUE:
        {
            // Fixed value: Uf = Ub
            return bc->fixedVectorValue();
        }

        case ZERO_GRADIENT:
        {
            return phi[face.ownerCell()];
        }

        case FIXED_GRADIENT:
        {
            Scalar dn = dot(face.dPf(), face.normal());

            if (bc->gradientType() == BCValueType::VECTOR)
            {
                return
                    phi[face.ownerCell()]
                  + bc->vectorGradient() * dn;
            }
            // Scalar gradient on a vector field
            Vector owner = phi[face.ownerCell()];
            Scalar grad = bc->scalarGradient();
            return
                Vector
                (
                    owner.x() + grad * dn,
                    owner.y() + grad * dn,
                    owner.z() + grad * dn
                );
        }

        default:
            throw
                std::runtime_error
                (
                    "Unknown BC type for face "
                  + std::to_string(face.idx())
                  + " in patch " + patch->patchName()
                  + ": "
                  + std::to_string(static_cast<int>(bc->type()))
                );
    }
}

void BoundaryConditions::linkFaces(std::vector<Face>& faces) const
{
    for (const auto& patch : patches_)
    {
        for
        (
            size_t faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            faces[faceIdx].setPatch(&patch);
        }
    }
}

std::string BoundaryConditions::bcTypeToString(BCType bctype)
{
    using enum BCType;
    switch (bctype)
    {
        case UNDEFINED: return "UNDEFINED";
        case FIXED_VALUE: return "FIXED_VALUE";
        case FIXED_GRADIENT: return "FIXED_GRADIENT";
        case ZERO_GRADIENT: return "ZERO_GRADIENT";
        case NO_SLIP: return "NO_SLIP";
        case K_WALL_FUNCTION: return "K_WALL_FUNCTION";
        case OMEGA_WALL_FUNCTION: return "OMEGA_WALL_FUNCTION";
        case NUT_WALL_FUNCTION: return "NUT_WALL_FUNCTION";
        default:
            throw
                std::runtime_error
                (
                    "Unknown BC type: "
                  + std::to_string(static_cast<int>(bctype))
                );
    }
}

void BoundaryConditions::validatePatchNames() const
{
    std::set<std::string> validNames;

    for (const auto& patch : patches_)
    {
        validNames.insert(patch.patchName());
    }

    for (const auto& entry : patchBoundaryData_)
    {
        if (!validNames.contains(entry.first))
        {
            std::string validList;
            for (const auto& name : validNames)
            {
                if (!validList.empty())
                {
                    validList += ", ";
                }
                validList += "'" + name + "'";
            }

            throw
                std::runtime_error
                (
                    "Boundary condition patch '"
                  + entry.first
                  + "' does not match any mesh patch. "
                    "Valid patch names: " + validList
                );
        }
    }
}

void BoundaryConditions::printSummary() const
{
    std::cout
        << std::endl
        << "--- Boundary Conditions Setup Summary ---" << std::endl;

    if (patches_.empty())
    {
        std::cout
            << "  No mesh patches loaded." << std::endl;

        return;
    }

    std::cout
        << "Total Mesh Patches Loaded: " << patches_.size()
        << std::endl;

    for (const auto& meshPatch : patches_)
    {
        std::cout
            << "  ------------------------------------"
            << std::endl;

        std::cout
            << "  Mesh Patch Name         : "
            << meshPatch.patchName() << std::endl;

        std::cout
            << "  Fluent Type             : "
            << meshPatch.fluentType() << std::endl;

        std::cout
            << "  Zone ID                 : "
            << meshPatch.zoneIdx() << std::endl;

        std::cout
            << "  Number of Faces         : "
            << meshPatch.numberOfBoundaryFaces() << std::endl;

        auto patchBCIterator =
            patchBoundaryData_.find(meshPatch.patchName());

        if
        (
            patchBCIterator != patchBoundaryData_.end()
         && !patchBCIterator->second.empty()
        )
        {
            std::cout
                << "  Configured Physical BCs :" << std::endl;

            for (const auto& fieldBCPair : patchBCIterator->second)
            {
                const std::string& fieldName = fieldBCPair.first;
                const BoundaryData& fbc = fieldBCPair.second;

                std::cout
                    << "      Field '" << fieldName
                    << "': Type: "
                    << bcTypeToString(fbc.type());

                if
                (
                    fbc.type() == BCType::FIXED_VALUE
                 || fbc.type() == BCType::NO_SLIP
                )
                {
                    std::cout
                        << ", Value: ";

                    if (fbc.valueType() == BCValueType::SCALAR)
                    {
                        std::cout
                            << fbc.scalarValue();
                    }
                    else if
                    (fbc.valueType() == BCValueType::VECTOR)
                    {
                        std::cout
                            << fbc.vectorValue();
                    }
                    else
                    {
                        throw
                            std::runtime_error
                            (
                                "Unknown BC value type"
                            );
                    }
                }
                else if (fbc.type() == BCType::FIXED_GRADIENT)
                {
                    std::cout
                        << ", Gradient: ";

                    if (fbc.gradientType() == BCValueType::SCALAR)
                    {
                        std::cout
                            << fbc.scalarGradient();
                    }
                    else if
                    (fbc.gradientType() == BCValueType::VECTOR)
                    {
                        std::cout
                            << fbc.vectorGradient();
                    }
                    else
                    {
                        throw
                            std::runtime_error
                            (
                                "Unknown BC gradient type: "
                              + std::to_string
                                (
                                    static_cast<int>
                                    (
                                        fbc.gradientType()
                                    )
                                )
                            );
                    }
                }
                else if (fbc.type() == BCType::ZERO_GRADIENT)
                {
                    std::cout
                        << " (implies zero gradient)";
                }
                else if
                (
                    fbc.type() == BCType::K_WALL_FUNCTION
                 || fbc.type() == BCType::OMEGA_WALL_FUNCTION
                 || fbc.type() == BCType::NUT_WALL_FUNCTION
                )
                {
                    std::cout
                        << " (wall function)";
                }

                std::cout
                    << std::endl;
            }
        }
    }

    std::cout
        << "  ------------------------------------" << std::endl;
}
