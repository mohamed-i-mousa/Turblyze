/******************************************************************************
 * @file BoundaryPatch.hpp
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
    VELOCITY_INLET,         ///< Velocity inlet boundary
    PRESSURE_INLET,         ///< Pressure inlet boundary
    PRESSURE_OUTLET,        ///< Pressure outlet boundary
    WALL,                   ///< Wall boundary
    SYMMETRY,               ///< Symmetry boundary
    PERIODIC,               ///< Periodic boundary
    MASS_FLOW_INLET,        ///< Mass flow inlet boundary
    OUTFLOW,                ///< Outflow boundary
    INTERFACE,              ///< Interface boundary
    INTERIOR,               ///< Interior boundary
    SOLID,                  ///< Solid boundary
    FLUID,                  ///< Fluid boundary
    UNDEFINED               ///< Undefined boundary type
};


class BoundaryPatch
{
public:

    /**
     * @brief Constructor for boundary patch
     * @param idx Zone identifier
     * @param startIdx Index of first face
     * @param endIdx Index of last face
     */
    BoundaryPatch
    (
        size_t idx,
        size_t startIdx,
        size_t endIdx
    ) noexcept
        : zoneIdx_(idx),
          firstFaceIdx_(startIdx),
          lastFaceIdx_(endIdx) {}

// Setter methods

    /**
     * @brief Set patch name
     * @param name New human-readable name
     */
    void setPatchName(std::string name) noexcept
    {
        patchName_ = std::move(name);
    }

    /**
     * @brief Set patch type
     * @param patchType New patch type
     */
    void setType(PatchType patchType) noexcept { type_ = patchType; }

// Accessor methods

    /**
     * @brief Get number of faces in this boundary patch
     * @return Number of boundary faces
     */
    size_t numberOfBoundaryFaces() const noexcept
    {
        return lastFaceIdx_ - firstFaceIdx_ + 1;
    }

    /**
     * @brief Get patch name
     * @return Patch name
     */
    const std::string& patchName() const noexcept { return patchName_; }

    /**
     * @brief Get patch type
     * @return Mapped patch type
     */
    PatchType type() const noexcept { return type_; }

    /**
     * @brief Get zone identifier
     * @return Zone ID from mesh file
     */
    size_t zoneIdx() const noexcept { return zoneIdx_; }

    /**
     * @brief Get first face index
     * @return Index of first face in patch
     */
    size_t firstFaceIdx() const noexcept { return firstFaceIdx_; }

    /**
     * @brief Get last face index
     * @return Index of last face in patch
     */
    size_t lastFaceIdx() const noexcept { return lastFaceIdx_; }

private:

// Private members

    /// Human-readable patch name
    std::string patchName_;

    /// Mapped boundary condition type
    PatchType type_ = PatchType::UNDEFINED;

    /// Zone identifier from mesh file
    size_t zoneIdx_;

    /// Index of first face in this patch
    size_t firstFaceIdx_;

    /// Index of last face in this patch
    size_t lastFaceIdx_;
};
