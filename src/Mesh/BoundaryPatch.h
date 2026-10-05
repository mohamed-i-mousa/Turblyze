/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryPatch.h
 * @brief Boundary patch representation and mesh connectivity management
 *
 * @details The BoundaryPatch class represents a set of faces on the domain
 * boundary. A patch is identified by a name (e.g., "inlet", "wall") and a
 * geometric zone ID from the mesh file.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <utility>

// Project headers
#include "Integer.h"
#include "StringTypes.h"

// *************************** enum class PatchType ***************************

enum class PatchType
{
    physical,               ///< Physical boundary patch
    processor               ///< Inter-rank cut in a decomposed mesh
};

// **************************** class BoundaryPatch ***************************

class BoundaryPatch
{
public:

// ************************* Special Member Functions *************************

    /// Constructor for boundary patch
    BoundaryPatch
    (
        Index idx,
        Index startIdx,
        Index endIdx
    ) noexcept
    :
        zoneIdx_(idx),
        firstFaceIdx_(startIdx),
        lastFaceIdx_(endIdx)
    {}

// ****************************** Setter Methods ******************************

    /// Set patch name
    void setName(Name patchName) noexcept
    {
        name_ = std::move(patchName);
    }

    /// Set patch type
    void setType(PatchType patchType) noexcept { type_ = patchType; }

// ***************************** Accessor Methods *****************************

    /// Get number of faces in this boundary patch
    [[nodiscard]] Count numFaces() const noexcept
    {
        return lastFaceIdx_ - firstFaceIdx_ + 1;
    }

    /// Get patch name
    [[nodiscard]] const Name& name() const noexcept
    {
        return name_;
    }

    /// Get patch type
    [[nodiscard]] PatchType type() const noexcept { return type_; }

    /// Get zone identifier
    [[nodiscard]] Index zoneIdx() const noexcept { return zoneIdx_; }

    /// Get first face index
    [[nodiscard]] Index firstFaceIdx() const noexcept
    {
        return firstFaceIdx_;
    }

    /// Get last face index
    [[nodiscard]] Index lastFaceIdx() const noexcept
    {
        return lastFaceIdx_;
    }

// ****************************** Private Members *****************************

private:

    /// Human-readable patch name
    Name name_;

    /// Patch classification: physical boundary or MPI processor cut
    PatchType type_ = PatchType::physical;

    /// Zone identifier from mesh file
    Index zoneIdx_;

    /// Index of first face in this patch
    Index firstFaceIdx_;

    /// Index of last face in this patch
    Index lastFaceIdx_;
};
