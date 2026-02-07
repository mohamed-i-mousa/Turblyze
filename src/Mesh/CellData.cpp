/******************************************************************************
 * @file CellData.cpp
 * @brief Implementation of cell-centered data containers
 *****************************************************************************/

#include "CellData.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>

// ****************************** Constructor *******************************

template<typename T>
CellData<T>::CellData
(
    const std::string& fieldName,
    size_t numCells
) : name_(fieldName),
    internalField_(numCells) {}

template<typename T>
CellData<T>::CellData
(
    const std::string& fieldName,
    size_t numCells,
    const T& initialValue
) : name_(fieldName),
    internalField_(numCells, initialValue) {}

template<typename T>
T& CellData<T>::operator[](size_t cellIndex)
{
    if (cellIndex >= internalField_.size())
    {
        throw   std::out_of_range
                (
                    "Cell index out of range in CellData '" + name_ + "'"
                );
    }

    return internalField_[cellIndex];
}

template<typename T>
const T& CellData<T>::operator[](size_t cellIndex) const
{
    if (cellIndex >= internalField_.size())
    {
        throw   std::out_of_range
                (
                    "Cell index out of range in CellData '" + name_ + "'"
                );
    }

    return internalField_[cellIndex];
}

// ***************************** Setter Methods *****************************

template<typename T>
void CellData<T>::setAll(const T& value)
{
    for (size_t i = 0; i < internalField_.size(); i++)
    {
        internalField_[i] = value;
    }
}

// ***************************** Utility Methods ****************************

template<typename T>
void CellData<T>::printSummary(size_t itemsToShow) const
{
    std::cout
        << "CellData: " << name_ << " (Size: " << internalField_.size() 
        << ")" << std::endl;

    for (size_t i = 0; i < std::min(internalField_.size(), itemsToShow); ++i)
    {
        std::cout
            << "  Cell " << i << ": " << internalField_[i] << std::endl;
    }

    if (internalField_.size() > itemsToShow)
    {
        std::cout
            << "  ..." << std::endl;
    }
}

// ************************* Template Instantiations *************************

template class CellData<Vector>;
template class CellData<Scalar>;