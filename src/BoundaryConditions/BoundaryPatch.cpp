/******************************************************************************
 * @file BoundaryPatch.cpp
 * @brief Implementation of boundary patch management
 *****************************************************************************/

#include "BoundaryPatch.hpp"
#include <iostream>


// ****************************** Free Functions ******************************

BoundaryConditionType mapFluentBCToEnum(const std::string& fluentType) 
{
    if (fluentType == "velocity-inlet") 
        return BoundaryConditionType::VELOCITY_INLET;

    if (fluentType == "pressure-inlet") 
        return BoundaryConditionType::PRESSURE_INLET;

    if (fluentType == "pressure-outlet") 
        return BoundaryConditionType::PRESSURE_OUTLET;

    if (fluentType == "wall") 
        return BoundaryConditionType::WALL;

    if (fluentType == "symmetry") 
        return BoundaryConditionType::SYMMETRY;

    if (fluentType == "periodic" || fluentType == "periodic-shadow") 
        return BoundaryConditionType::PERIODIC;

    if (fluentType == "mass-flow-inlet") 
        return BoundaryConditionType::MASS_FLOW_INLET;

    if (fluentType == "outflow") 
        return BoundaryConditionType::OUTFLOW;

    if (fluentType == "interface") 
        return BoundaryConditionType::INTERFACE;

    if (fluentType == "interior") 
        return BoundaryConditionType::INTERIOR;

    if (fluentType == "solid") 
        return BoundaryConditionType::SOLID;

    if (fluentType == "fluid") 
        return BoundaryConditionType::FLUID;

    std::cerr
        << "Warning: Unknown Fluent boundary type encountered: "
        << fluentType << std::endl;
                
    return BoundaryConditionType::UNDEFINED;
}


// ****************************** Public Methods ******************************

size_t BoundaryPatch::numberOfBoundaryFaces() const 
{
    return lastFaceIdx_ - firstFaceIdx_ + 1;
}