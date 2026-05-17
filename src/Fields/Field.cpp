/******************************************************************************
 * @file Field.cpp
 * @brief Implementation of Field identifier helpers
 *****************************************************************************/

#include "Field.h"

#include "ErrorHandler.h"


std::string_view fieldToString(Field field)
{
    using enum Field;
    switch (field)
    {
        case Ux:    return "Ux";
        case Uy:    return "Uy";
        case Uz:    return "Uz";
        case p:     return "p";
        case pCorr: return "pCorr";
        case k:     return "k";
        case omega: return "omega";
        case nut:   return "nut";
    }

    FatalError("Corrupted Field value");
}
