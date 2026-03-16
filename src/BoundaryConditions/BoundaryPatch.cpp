/******************************************************************************
 * @file BoundaryPatch.cpp
 * @brief Implementation of boundary patch management
 *****************************************************************************/

#include "BoundaryPatch.hpp"

#include <iostream>


// *************************** Static Methods *********************************

BoundaryConditionType BoundaryPatch::mapFluentBCToEnum
(
    std::string_view fluentType
)
{
    for (const auto& [name, type] : bcMappings_)
    {
        if (fluentType == name)
            return type;
    }

    std::cerr
        << "Warning: Unknown Fluent boundary type encountered: "
        << fluentType << std::endl;

    return BoundaryConditionType::UNDEFINED;
}
