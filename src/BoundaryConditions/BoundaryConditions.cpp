/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryConditions.cpp
 * @brief Implementation of boundary conditions management system
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "BoundaryConditions.h"

// Standard library headers
#include <iostream>
#include <utility>
#include <set>

// Project headers
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
    const Name& patchName,
    Field field,
    Scalar value
)
{
    BoundaryData bcData;
    bcData.setFixedValue(value);
    setBC(patchName, field, std::move(bcData));
}


void BoundaryConditions::setFixedGradient
(
    const Name& patchName,
    Field field,
    Scalar gradient
)
{
    BoundaryData bcData;
    bcData.setFixedGradient(gradient);
    setBC(patchName, field, std::move(bcData));
}


void BoundaryConditions::setZeroGradient
(
    const Name& patchName,
    Field field
)
{
    BoundaryData bcData;
    bcData.setZeroGradient();
    setBC(patchName, field, std::move(bcData));
}


void BoundaryConditions::setNoSlip
(
    const Name& patchName,
    Field field
)
{
    BoundaryData bcData;
    bcData.setNoSlip();
    setBC(patchName, field, std::move(bcData));
}


void BoundaryConditions::setWallFunctionType
(
    const Name& patchName,
    Field field,
    BCType wallType
)
{
    BoundaryData bcData;
    bcData.setWallFunctionType(wallType);
    setBC(patchName, field, std::move(bcData));
}

// ***************************** Accessor Methods *****************************

const BoundaryPatch& BoundaryConditions::patch(const Name& patchName) const
{
    for (const auto& patch : patches_)
    {
        if (patch.patchName() == patchName)
        {
            return patch;
        }
    }

    FatalError("Patch " + patchName + " not found");
}


const BoundaryData& BoundaryConditions::fieldBC
(
    const Name& patchName,
    Field field
) const
{
    const auto patchIterator = patchBoundaryData_.find(patchName);

    if (patchIterator != patchBoundaryData_.end())
    {
        const auto& fieldMap = patchIterator->second;

        const auto fieldIterator = fieldMap.find(field);

        if (fieldIterator != fieldMap.end())
        {
            return fieldIterator->second;
        }
    }

    FatalError
    (
        "Boundary condition not found for patch " + patchName
        + " and field " + Name(fieldToString(field))
    );
}


Scalar BoundaryConditions::boundaryFaceValue
(
    Field field,
    const ScalarField& phi,
    const Face& boundaryFace
) const
{
    if (!boundaryFace.patch().has_value())
    {
        FatalError
        (
            "Face " + std::to_string(boundaryFace.idx())
            + " has no linked patch. "
              "Ensure linkFaces() is called before solving."
        );
    }

    const BoundaryPatch& patch = boundaryFace.patch()->get();
    const BoundaryData& bc = fieldBC(patch.patchName(), field);

    using enum BCType;
    switch (bc.type())
    {
        case noSlip:
        {
            return S(0.0);
        }
        
        case fixedValue:
        {
            // Fixed value: φf = φb
            return bc.fixedScalarValue();
        }

        case kWallFunction:
        case omegaWallFunction:
        case nutWallFunction:
        case zeroGradient:
        {
            // Zero gradient: φf = φP
            return phi[boundaryFace.ownerCell()];
        }

        case fixedGradient:
        {
            // Fixed gradient: φf = φP + grad * distance
            const Scalar dn = dot(boundaryFace.dPf(), boundaryFace.normal());
            
            return
                phi[boundaryFace.ownerCell()]
              + bc.fixedScalarGradient() * dn;
        }

        case undefined:
        default:
        {
            FatalError
            (
                "Corrupted BCType value for face "
              + std::to_string(boundaryFace.idx())
              + " in patch " + patch.patchName()
            );
        }

    }
}


void BoundaryConditions::linkFaces(MutableFaceListRef faces)
{
    for (const auto& patch : patches_)
    {
        if (patch.firstFaceIdx() > patch.lastFaceIdx())
        {
            FatalError
            (
                "Boundary patch '" + patch.patchName()
              + "' has an inverted face range ["
              + std::to_string(patch.firstFaceIdx()) + ", "
              + std::to_string(patch.lastFaceIdx()) + "]."
            );
        }

        if (patch.lastFaceIdx() >= faces.size())
        {
            FatalError
            (
                "Boundary patch '" + patch.patchName()
              + "' references face index "
              + std::to_string(patch.lastFaceIdx())
              + " outside the valid range [0, "
              + std::to_string(faces.size()) + ")."
            );
        }

        for
        (
            Index faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            faces[faceIdx].setPatch(patch);
        }
    }

    linked_ = true;
}


void BoundaryConditions::validatePatchNames() const
{
    // std::set guarantees uniqueness
    std::set<Name> validNames;

    for (const auto& patch : patches_)
    {
        validNames.insert(patch.patchName());
    }

    for (const auto& entry : patchBoundaryData_)
    {
        if (!validNames.contains(entry.first))
        {
            Message validList;
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

        const auto patchBCIterator =
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
                const Field field = fieldBCPair.first;
                const BoundaryData& fbc = fieldBCPair.second;

                std::cout
                    << "      Field '" << fieldToString(field)
                    << "': Type: "
                    << bcTypeToString(fbc.type());

                using enum BCType;
                switch (fbc.type())
                {
                    case fixedValue:
                    case noSlip:
                        std::cout
                            << ", Value: "
                            << fbc.fixedScalarValue();
                        break;

                    case fixedGradient:
                        std::cout
                            << ", Gradient: "
                            << fbc.fixedScalarGradient();
                        break;

                    case zeroGradient:
                        std::cout
                            << " (implies zero gradient)";
                        break;

                    case kWallFunction:
                    case omegaWallFunction:
                    case nutWallFunction:
                        std::cout
                            << " (wall function)";
                        break;

                    case undefined:
                        break;
                }

                std::cout
                    << '\n';
            }
        }
    }

    std::cout
        << "  ------------------------------------" << '\n';
}
