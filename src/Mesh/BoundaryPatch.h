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

#include <string>


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


class BoundaryPatch
{
public:

    /// Constructor for boundary patch
    BoundaryPatch
    (
        size_t idx,
        size_t startIdx,
        size_t endIdx
    ) noexcept
    :
        zoneIdx_(idx),
        firstFaceIdx_(startIdx),
        lastFaceIdx_(endIdx)
    {}

// Setter methods

    /// Set patch name
    void setPatchName(std::string name) noexcept
    {
        patchName_ = std::move(name);
    }

    /// Set patch type
    void setType(PatchType patchType) noexcept { type_ = patchType; }

// Accessor methods

    /// Get number of faces in this boundary patch
    [[nodiscard]] size_t numBoundaryFaces() const noexcept
    {
        return lastFaceIdx_ - firstFaceIdx_ + 1;
    }

    /// Get patch name
    [[nodiscard]] const std::string& patchName() const noexcept
    {
        return patchName_;
    }

    /// Get patch type
    [[nodiscard]] PatchType type() const noexcept { return type_; }

    /// Get zone identifier
    [[nodiscard]] size_t zoneIdx() const noexcept { return zoneIdx_; }

    /// Get first face index
    [[nodiscard]] size_t firstFaceIdx() const noexcept
    {
        return firstFaceIdx_;
    }

    /// Get last face index
    [[nodiscard]] size_t lastFaceIdx() const noexcept
    {
        return lastFaceIdx_;
    }

private:

// Private members

    /// Human-readable patch name
    std::string patchName_;

    /// Mapped boundary condition type
    PatchType type_ = PatchType::undefined;

    /// Zone identifier from mesh file
    size_t zoneIdx_;

    /// Index of first face in this patch
    size_t firstFaceIdx_;

    /// Index of last face in this patch
    size_t lastFaceIdx_;
};
