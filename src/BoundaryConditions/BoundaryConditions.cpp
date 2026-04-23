/******************************************************************************
 * @file BoundaryConditions.cpp
 * @brief Implementation of boundary conditions management system
 *****************************************************************************/

#include "BoundaryConditions.h"

#include <iostream>
#include <utility>
#include <set>

#include "ErrorHandler.h"


// ****************************** Setter Methods ******************************

void BoundaryConditions::addPatch(BoundaryPatch patch)
{
    if (linked_)
    {
        FatalError
        (
            "Cannot add patch after linkFaces() has been called — "
            "stored face pointers would become invalid."
        );
    }
    patches_.push_back(std::move(patch));
}

void BoundaryConditions::setFixedValue
(
    const std::string& patchName,
    const std::string& fieldName,
    Scalar value
)
{
    BoundaryData bcData;
    bcData.setFixedValue(value);
    setBC(patchName, fieldName, std::move(bcData));
}

void BoundaryConditions::setFixedValue
(
    const std::string& patchName,
    const std::string& fieldName,
    const Vector& value
)
{
    BoundaryData bcData;
    bcData.setFixedValue(value);
    setBC(patchName, fieldName, std::move(bcData));
}

void BoundaryConditions::setFixedGradient
(
    const std::string& patchName,
    const std::string& fieldName,
    Scalar gradient
)
{
    BoundaryData bcData;
    bcData.setFixedGradient(gradient);
    setBC(patchName, fieldName, std::move(bcData));
}

void BoundaryConditions::setZeroGradient
(
    const std::string& patchName,
    const std::string& fieldName
)
{
    BoundaryData bcData;
    bcData.setZeroGradient();
    setBC(patchName, fieldName, std::move(bcData));
}

void BoundaryConditions::setNoSlip
(
    const std::string& patchName,
    const std::string& fieldName
)
{
    BoundaryData bcData;
    bcData.setNoSlip();
    setBC(patchName, fieldName, std::move(bcData));
}

void BoundaryConditions::setWallFunctionType
(
    const std::string& patchName,
    const std::string& fieldName,
    BCType wallType
)
{
    BoundaryData bcData;
    bcData.setWallFunctionType(wallType);
    setBC(patchName, fieldName, std::move(bcData));
}


// ****************************** Accessor Methods ******************************

const BoundaryPatch& BoundaryConditions::patch(const std::string& name) const
{
    for (const auto& patch : patches_)
    {
        if (patch.patchName() == name)
        {
            return patch;
        }
    }

    FatalError("Patch " + name + " not found");
}

const BoundaryData& BoundaryConditions::fieldBC
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
            return fieldIterator->second;
        }
    }

    FatalError
    (
        "Boundary condition not found for patch " + patchName
        + " and field " + fieldName
    );
}

Scalar BoundaryConditions::boundaryFaceValue
(
    const std::string& fieldName,
    const ScalarField& phi,
    const Face& face,
    std::optional<int> componentIdx
) const
{
    if (!face.patch().has_value())
    {
        FatalError
        (
            "Face " + std::to_string(face.idx())
            + " has no linked patch. "
              "Ensure linkFaces() is called before solving."
        );
    }

    const BoundaryPatch& patch = face.patch()->get();
    const BoundaryData& bc = fieldBC(patch.patchName(), fieldName);

    using enum BCType;
    switch (bc.type())
    {
        case NO_SLIP:
        {
            return S(0.0);
        }
        
        case FIXED_VALUE:
        {
            // Extract scalar from vector BC using component index
            if (componentIdx && bc.valueType() == BCValueType::VECTOR)
            {
                const Vector& v = bc.fixedVectorValue();
                switch (*componentIdx)
                {
                    case 0: return v.x();
                    case 1: return v.y();
                    case 2: return v.z();
                    default:
                        FatalError
                        (
                            "Invalid component index "
                            + std::to_string(*componentIdx)
                            + " for vector BC on face "
                            + std::to_string(face.idx())
                        );
                }
            }

            // Fixed value: φf = φb
            return bc.fixedScalarValue();
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
              + bc.fixedScalarGradient() * dn;
        }

        case UNDEFINED:
        default:
        {
            FatalError
            (
                "Corrupted BCType value for face "
              + std::to_string(face.idx())
              + " in patch " + patch.patchName()
            );
        }

    }
}

Vector BoundaryConditions::boundaryVectorFaceValue
(
    const std::string& fieldName,
    const VectorField& phi,
    const Face& face
) const
{
    if (!face.patch().has_value())
    {
        FatalError
        (
            "Face " + std::to_string(face.idx())
            + " has no linked patch. "
              "Ensure linkFaces() is called before solving."
        );
    }

    const BoundaryPatch& patch = face.patch()->get();
    const BoundaryData& bc = fieldBC(patch.patchName(), fieldName);

    using enum BCType;
    switch (bc.type())
    {
        case NO_SLIP:
        {
            return Vector{};
        }

        case FIXED_VALUE:
        {
            // Fixed value: Uf = Ub
            return bc.fixedVectorValue();
        }

        case K_WALL_FUNCTION:
        case OMEGA_WALL_FUNCTION:
        case NUT_WALL_FUNCTION:
        case ZERO_GRADIENT:
        {
            return phi[face.ownerCell()];
        }

        case FIXED_GRADIENT:
        {
            // Fixed gradient: Uf = UP + grad * distance (per component)
            Scalar dn = dot(face.dPf(), face.normal());
            Scalar grad = bc.scalarGradient();
            Vector owner = phi[face.ownerCell()];
            return
                Vector
                (
                    owner.x() + grad * dn,
                    owner.y() + grad * dn,
                    owner.z() + grad * dn
                );
        }

        case UNDEFINED:
        default:
        {
            FatalError
            (
                "Corrupted BCType value for face "
              + std::to_string(face.idx())
              + " in patch " + patch.patchName()
            );
        }
    }
}

void BoundaryConditions::linkFaces(std::span<Face> faces)
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
            faces[faceIdx].setPatch(patch);
        }
    }
    
    linked_ = true;
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
    }

    FatalError("Corrupted BCType value");
}

void BoundaryConditions::validatePatchNames() const
{
    // std::set guarantees uniqueness
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

            FatalError
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
        << '\n'
        << "--- Boundary Conditions Setup Summary ---" << '\n';

    if (patches_.empty())
    {
        std::cout
            << "  No mesh patches loaded." << '\n';

        return;
    }

    std::cout
        << "Total Mesh Patches Loaded: " << patches_.size()
        << '\n';

    for (const auto& meshPatch : patches_)
    {
        std::cout
            << "  ------------------------------------"
            << '\n';

        std::cout
            << "  Mesh Patch Name         : "
            << meshPatch.patchName() << '\n';

        std::cout
            << "  Zone ID                 : "
            << meshPatch.zoneIdx() << '\n';

        std::cout
            << "  Number of Faces         : "
            << meshPatch.numBoundaryFaces() << '\n';

        auto patchBCIterator =
            patchBoundaryData_.find(meshPatch.patchName());

        if
        (
            patchBCIterator != patchBoundaryData_.end()
         && !patchBCIterator->second.empty()
        )
        {
            std::cout
                << "  Configured Physical BCs :" << '\n';

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
                    else if (fbc.valueType() == BCValueType::VECTOR)
                    {
                        std::cout
                            << fbc.vectorValue();
                    }
                    else
                    {
                        FatalError("Unknown BC value type");
                    }
                }
                else if (fbc.type() == BCType::FIXED_GRADIENT)
                {
                    std::cout
                        << ", Gradient: "
                        << fbc.scalarGradient();
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
                    << '\n';
            }
        }
    }

    std::cout
        << "  ------------------------------------" << '\n';
}
