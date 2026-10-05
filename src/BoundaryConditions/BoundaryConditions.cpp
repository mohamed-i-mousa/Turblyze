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
#include <format>
#include <iostream>

// Project headers
#include "BoundaryPatch.h"
#include "ErrorHandler.h"
#include "Face.h"
#include "Field.h"
#include "Mesh.h"

// ************************* Special Member Functions *************************

BoundaryConditions::BoundaryConditions
(
    BCs boundaryConditions,
    const Mesh& mesh
)
:
    boundaryConditions_(std::move(boundaryConditions))
{
    for (const BoundaryPatch& patch : mesh.patches())
    {
        if (patch.type() == PatchType::processor)
        {
            continue;
        }

        const auto patchIt = boundaryConditions_.find(patch.name());
        if (patchIt == boundaryConditions_.end())
        {
            continue;
        }

        for (auto& [field, bcPtr] : patchIt->second)
        {
            if (bcPtr != nullptr)
            {
                bcPtr->updateCoeffs(mesh, patch);
            }
        }
    }
}

// ***************************** Accessor Methods *****************************

const BoundaryType& BoundaryConditions::boundaryType
(
    const Name& patchName,
    Field field
) const
{
    const BoundaryType* bc = find(patchName, field);

    if (bc == nullptr)
    {
        missing(patchName, field);
    }

    return *bc;
}


const BoundaryType& BoundaryConditions::boundaryType
(
    const BoundaryPatch& patch,
    Field field
) const
{
    return boundaryType(patch.name(), field);
}


const BoundaryType& BoundaryConditions::boundaryType
(
    const Face& face,
    Field field
) const
{
    const BoundaryPatch* patch = face.patch();
    if (patch == nullptr || patch->type() == PatchType::processor)
    {
        FatalError("Face is not on a physical boundary patch");
    }

    return boundaryType(patch->name(), field);
}


Scalar BoundaryConditions::faceValue
(
    const Face& face,
    Scalar ownerValue,
    Field field
) const noexcept
{
    const BoundaryPatch* patch = face.patch();
    if (patch == nullptr || patch->type() == PatchType::processor)
    {
        return ownerValue;
    }

    const BoundaryType* bc = find(patch->name(), field);
    if (bc == nullptr)
    {
        return ownerValue;
    }

    const Index localIdx = face.idx() - patch->firstFaceIdx();
    return bc->faceValue(localIdx, ownerValue);
}

// ****************************** Public Methods ******************************

void BoundaryConditions::refresh
(
    const Mesh& mesh,
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz
)
{
    for (const BoundaryPatch& patch : mesh.patches())
    {
        if (patch.type() == PatchType::processor)
        {
            continue;
        }

        const auto patchIt = boundaryConditions_.find(patch.name());
        if (patchIt == boundaryConditions_.end())
        {
            missing(patch.name(), Field::Ux);
        }

        for (const Field field : {Field::Ux, Field::Uy, Field::Uz})
        {
            const auto fieldIt = patchIt->second.find(field);
            if (fieldIt == patchIt->second.end() || fieldIt->second == nullptr)
            {
                missing(patch.name(), field);
            }

            fieldIt->second->refreshCoeffs(mesh, patch, Ux, Uy, Uz);
        }
    }
}


void BoundaryConditions::printSummary() const
{
    std::cout << "\n--- Boundary Conditions Setup Summary ---\n";

    if (boundaryConditions_.empty())
    {
        std::cout << "  No boundary conditions loaded.\n";
        return;
    }

    std::cout << std::format
    (
        "Total Patches Configured: {}\n",
        boundaryConditions_.size()
    );

    for (const auto& [patchName, fieldMap] : boundaryConditions_)
    {
        std::cout << std::format
        (
            "  ------------------------------------\n"
            "  Patch Name              : {}\n"
            "  Configured Fields       :\n",
            patchName
        );

        for (const auto& [field, bc] : fieldMap)
        {
            std::cout << std::format
            (
                "      Field '{}': Type: ",
                fieldToString(field)
            );

            if (bc != nullptr)
            {
                std::cout << bc->typeName();
            }

            std::cout << '\n';
        }
    }

    std::cout << "  ------------------------------------\n";
}

// ****************************** Private Methods *****************************

const BoundaryType* BoundaryConditions::find
(
    const Name& patchName,
    Field field
) const noexcept
{
    const auto patchIt = boundaryConditions_.find(patchName);
    if (patchIt == boundaryConditions_.end())
    {
        return nullptr;
    }

    const auto fieldIt = patchIt->second.find(field);
    if (fieldIt == patchIt->second.end())
    {
        return nullptr;
    }

    return fieldIt->second.get();
}


void BoundaryConditions::missing(const Name& patchName, Field field)
{
    FatalError
    (
        "Boundary condition not found for patch '" + patchName
      + "' and field '" + Name(fieldToString(field)) + "'."
    );
}
