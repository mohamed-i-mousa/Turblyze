/******************************************************************************
 * @file BoundaryPatch.h
 * @brief Boundary patch representation and mesh connectivity management
 * 
 * @class BoundaryPatch
 * 
 * Represents a collection of boundary faces sharing the same boundary
 * condition type and patch name. Manages mesh connectivity information,
 * face-to-patch mapping, and boundary type classification for CFD solver.
 * 
 * Key features:
 * - Face grouping by boundary condition type and mesh patch name
 * - Mesh connectivity storage (face IDs, node indices, element types)
 * - Boundary type classification from mesh file input
 * - Integration with BoundaryConditions management system
 *****************************************************************************/

#ifndef BOUNDARYPATCH_H
#define BOUNDARYPATCH_H

#include <string>
#include <vector>
#include <map>
#include <stdexcept>

#include "Face.hpp"

/**
 * @enum BoundaryConditionType
 * @brief Enumeration of boundary condition types from mesh files
 */
enum class BoundaryConditionType 
{
    VELOCITY_INLET,     ///< Velocity inlet boundary
    PRESSURE_INLET,     ///< Pressure inlet boundary
    PRESSURE_OUTLET,    ///< Pressure outlet boundary
    WALL,               ///< Wall boundary
    SYMMETRY,           ///< Symmetry boundary
    PERIODIC,           ///< Periodic boundary
    MASS_FLOW_INLET,    ///< Mass flow inlet boundary
    OUTFLOW,            ///< Outflow boundary
    INTERFACE,          ///< Interface boundary
    INTERIOR,           ///< Interior boundary
    SOLID,              ///< Solid boundary
    FLUID,              ///< Fluid boundary
    UNDEFINED           ///< Undefined boundary type
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
     * @param id Zone identifier
     * @param start_id Index of first face
     * @param end_id Index of last face
     */
    BoundaryPatch
    (
        size_t id, 
        size_t start_id, 
        size_t end_id
    ) : zoneID_(id), 
        firstFaceIdx_(start_id), 
        lastFaceIdx_(end_id) {}

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
    size_t numberOfBoundaryFaces() const;

    /** 
     * @brief Get patch name 
     * @return Human-readable patch name 
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
    size_t zoneID() const { return zoneID_; }
    
    /** 
     * @brief Get first face index 
     * @return Index of first face in patch 
     */
    size_t firstFaceIndex() const { return firstFaceIdx_; }
    
    /** 
     * @brief Get last face index 
     * @return Index of last face in patch
     */
    size_t lastFaceIndex() const { return lastFaceIdx_; }
    
private:

// Private members

    /// Human-readable patch name
    std::string patchName_;
    
    /// Original Fluent boundary type string
    std::string fluentType_;
    
    /// Mapped boundary condition type
    BoundaryConditionType type_;
    
    /// Zone identifier from mesh file
    size_t zoneID_;
    
    /// Index of first face in this patch
    size_t firstFaceIdx_;
    
    /// Index of last face in this patch
    size_t lastFaceIdx_;
};

#endif