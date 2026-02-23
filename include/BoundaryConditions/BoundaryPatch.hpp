/******************************************************************************
 * @file BoundaryPatch.hpp
 * @brief Boundary patch representation and mesh connectivity management
 *
 * @details This header defines the BoundaryPatch class, which represents a 
 * set of faces on the domain boundary. A patch is identified by a name
 * (e.g., "inlet", "wall") and a geometric zone ID from the mesh file.
 * 
 * @class BoundaryPatch
 * 
 * The BoundaryPatch class provides:
 * - Identification of boundary zones (name, ID, type)
 * - Topological range definitions (start face index, end face index)
 * - Mapping between mesh file types (e.g., Fluent strings) and internal enums
 * - Helper methods for querying patch size and face validity
 *****************************************************************************/

#ifndef BOUNDARY_PATCH_HPP
#define BOUNDARY_PATCH_HPP

#include <string>
#include <vector>
#include <map>
#include <stdexcept>

#include "Face.hpp"

/**
 * @enum BoundaryConditionType
 * @brief Enumeration of boundary condition types
 */
enum class BoundaryConditionType
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

/**
 * @brief Maps Fluent boundary type string to enumeration
 * @param fluentType String representation from Fluent mesh file
 * @return Corresponding BoundaryConditionType enumeration
 */
BoundaryConditionType mapFluentBCToEnum(const std::string& fluentType);


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
    ) : zoneIdx_(idx), 
        firstFaceIdx_(startIdx), 
        lastFaceIdx_(endIdx) {}

// Setter methods
    
    /** 
     * @brief Set patch name
     * @param name New human-readable name 
     */
    void setPatchName(const std::string& name) { patchName_ = name; }
    
    /** 
     * @brief Set Fluent type 
     * @param typeStr New Fluent boundary type string 
     */
    void setFluentType(const std::string& typeStr) { fluentType_ = typeStr; }
    
    /** 
     * @brief Set boundary condition type 
     * @param bcType New boundary condition type 
     */
    void setType(BoundaryConditionType bcType) { type_ = bcType; }

// Accessor methods

    /**
     * @brief Get number of faces in this boundary patch
     * @return Number of boundary faces
     */
    size_t numberOfBoundaryFaces() const
    { 
        return lastFaceIdx_ - firstFaceIdx_ + 1; 
    }

    /** 
     * @brief Get patch name 
     * @return Patch name 
     */
    const std::string& patchName() const { return patchName_; }
    
    /** 
     * @brief Get Fluent type string 
     * @return Original Fluent boundary type 
     */
    const std::string& fluentType() const { return fluentType_; }
    
    /** 
     * @brief Get boundary condition type 
     * @return Mapped boundary condition type 
     */
    BoundaryConditionType type() const { return type_; }
    
    /** 
     * @brief Get zone identifier 
     * @return Zone ID from mesh file 
     */
    size_t zoneIdx() const { return zoneIdx_; }
    
    /** 
     * @brief Get first face index 
     * @return Index of first face in patch 
     */
    size_t firstFaceIdx() const { return firstFaceIdx_; }
    
    /** 
     * @brief Get last face index 
     * @return Index of last face in patch
     */
    size_t lastFaceIdx() const { return lastFaceIdx_; }

private:

// Private members

    /// Human-readable patch name
    std::string patchName_;
    
    /// Original Fluent boundary type string
    std::string fluentType_;
    
    /// Mapped boundary condition type
    BoundaryConditionType type_;
    
    /// Zone identifier from mesh file
    size_t zoneIdx_;
    
    /// Index of first face in this patch
    size_t firstFaceIdx_;
    
    /// Index of last face in this patch
    size_t lastFaceIdx_;
};

#endif // BOUNDARY_PATCH_HPP
