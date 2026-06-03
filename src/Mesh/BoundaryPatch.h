/******************************************************************************
 * @file BoundaryPatch.h
 * @brief Boundary patch representation and mesh connectivity management
 *
 * @details This header defines the BoundaryPatch class, which represents a
 * set of faces on the domain boundary. A patch is identified by a name
 * (e.g., "inlet", "wall") and a geometric zone ID from the mesh file.
 *
 * @enum PatchType
 * - Enumeration of boundary condition types
 *
 * @class BoundaryPatch
 * - Identification of boundary zones (name, ID, type)
 * - Topological range definitions (start face index, end face index)
 * - Helper methods for querying patch size and face validity
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <utility>

#include "Integer.h"
#include "StringTypes.h"

// *************************** enum class PatchType ***************************

enum class PatchType
{
    velocityInlet,          ///< Velocity inlet boundary
    pressureInlet,          ///< Pressure inlet boundary
    pressureOutlet,         ///< Pressure outlet boundary
    wall,                   ///< Wall boundary
    symmetry,               ///< Symmetry boundary
    periodic,               ///< Periodic boundary
    massFlowInlet,          ///< Mass flow inlet boundary
    outflow,                ///< Outflow boundary
    interface,              ///< Interface boundary
    interior,               ///< Interior boundary
    solid,                  ///< Solid boundary
    fluid,                  ///< Fluid boundary
    undefined               ///< Undefined boundary type
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
    void setPatchName(Name patchName) noexcept
    {
        patchName_ = std::move(patchName);
    }

    /// Set patch type
    void setType(PatchType patchType) noexcept { type_ = patchType; }

// ***************************** Accessor Methods *****************************

    /// Get number of faces in this boundary patch
    [[nodiscard]] Count numBoundaryFaces() const noexcept
    {
        return lastFaceIdx_ - firstFaceIdx_ + 1;
    }

    /// Get patch name
    [[nodiscard]] const Name& patchName() const noexcept
    {
        return patchName_;
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
    Name patchName_;

    /// Mapped boundary condition type
    PatchType type_ = PatchType::undefined;

    /// Zone identifier from mesh file
    Index zoneIdx_;

    /// Index of first face in this patch
    Index firstFaceIdx_;

    /// Index of last face in this patch
    Index lastFaceIdx_;
};
