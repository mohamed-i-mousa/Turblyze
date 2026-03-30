/******************************************************************************
 * @file CellData.cpp
 * @brief Implementation of cell-centered data containers
 *****************************************************************************/

#include "CellData.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>

// ****************************** Constructor *******************************

template<typename T>
CellData<T>::CellData
(
    std::string fieldName,
    size_t numCells
) : name_(std::move(fieldName)),
    internalField_(numCells)
{
    assert(numCells > 0);
}

template<typename T>
CellData<T>::CellData
(
    std::string fieldName,
    size_t numCells,
    const T& initialValue
) : name_(std::move(fieldName)),
    internalField_(numCells, initialValue)
{
    assert(numCells > 0);
}

// ***************************** Setter Methods *****************************

template<typename T>
void CellData<T>::setAll(const T& value)
{
    std::fill(internalField_.begin(), internalField_.end(), value);
}

// ***************************** Utility Methods ****************************

template<typename T>
void CellData<T>::printSummary(size_t itemsToShow) const
{
    std::cout
        << "CellData: " << name_ << " (Size: " << internalField_.size()
        << ")\n";

    for (size_t i = 0; i < std::min(internalField_.size(), itemsToShow); ++i)
    {
        std::cout
            << "  Cell " << i << ": " << internalField_[i] << '\n';
    }

    if (internalField_.size() > itemsToShow)
    {
        std::cout
            << "  ...\n";
    }
}

// ************************* Template Instantiations *************************

template class CellData<Vector>;
template class CellData<Scalar>;