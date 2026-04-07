/******************************************************************************
 * @file CellData.cpp
 * @brief Implementation of cell-centered data containers
 *****************************************************************************/

#include "CellData.hpp"

#include <iostream>

#include "ErrorHandler.hpp"


// ***************************** Helper Methods ******************************

template<typename T>
void CellData<T>::printSummary(size_t itemsToShow) const
{
    std::cout
        << "CellData (Size: " << internalField_.size() << ")\n";

    for
    (
        size_t cellIdx = 0;
        cellIdx < std::min(internalField_.size(),itemsToShow);
        ++cellIdx
    )
    {
        std::cout
            << "  Cell " << cellIdx << ": " << internalField_[cellIdx] << '\n';
    }

    if (internalField_.size() > itemsToShow)
    {
        std::cout
            << "  ...\n";
    }
}

// ************************* Template Instantiations **************************

template class CellData<Vector>;
template class CellData<Scalar>;
