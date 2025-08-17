#ifndef BOUNDARYPATCH_H
#define BOUNDARYPATCH_H

#include <string>
#include <vector>
#include <map>

#include "Face.h"

/**
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

/**
 * @brief Represents a boundary patch in the mesh
 * 
 * A boundary patch is a collection of boundary faces that share
 * the same boundary condition type. Contains mesh connectivity
 * information and boundary type classification.
 */
struct BoundaryPatch 
{
    /// Human-readable patch name
    std::string patchName;
    
    /// Original Fluent boundary type string
    std::string fluentType;
    
    /// Mapped boundary condition type
    BoundaryConditionType type;
    
    /// Zone identifier from mesh file
    size_t zoneID;
    
    /// Index of first face in this patch
    size_t firstFaceIndex;
    
    /// Index of last face in this patch
    size_t lastFaceIndex;

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
    ) : zoneID(id), 
        firstFaceIndex(start_id), 
        lastFaceIndex(end_id) {}

    /**
     * @brief Get number of faces in this boundary patch
     * @return Number of boundary faces
     */
    size_t getNumberOfBoundaryFaces() const;
};

#endif