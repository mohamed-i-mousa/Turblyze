#include "CellData.h"
#include <stdexcept>
#include <iostream>
#include <algorithm>

template<typename T>
CellData<T>::CellData
(
    const std::string& fieldName,
    size_t numCells
) : name(fieldName),
    internalField(numCells) {}

template<typename T>
CellData<T>::CellData
(
    const std::string& fieldName,
    size_t numCells,
    const T& initialValue
) : name(fieldName),
    internalField(numCells, initialValue) {}

template<typename T>
T& CellData<T>::operator[](size_t cellIndex)
{
    if (cellIndex >= internalField.size())
    {
        throw std::out_of_range
        (
            "Cell index out of range in CellData '" + name + "'"
        );
    }

    return internalField[cellIndex];
}

template<typename T>
const T& CellData<T>::operator[](size_t cellIndex) const
{
    if (cellIndex >= internalField.size())
    {
        throw std::out_of_range
        (
            "Cell index out of range in CellData '" + name + "'"
        );
    }

    return internalField[cellIndex];
}


template<typename T>
void CellData<T>::setAll(const T& value) 
{
    for (size_t i = 0; i < internalField.size(); i++)
    {
        internalField[i] = value;
    }
}

template<typename T>
void CellData<T>::printSummary(size_t itemsToShow) const
{
    std::cout   << "CellData: " << name << " (Size: " 
                << internalField.size() << ")" << std::endl;

    for (size_t i = 0; i < std::min(internalField.size(), itemsToShow); ++i)
    {
        std::cout   << "  Cell " << i << ": " << internalField[i] 
                    << std::endl;
    }
    
    if (internalField.size() > itemsToShow)
    {
        std::cout << "  ..." << std::endl;
    }
}

template class CellData<Vector>;
template class CellData<Scalar>;
