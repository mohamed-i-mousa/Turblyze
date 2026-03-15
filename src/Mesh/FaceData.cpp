/******************************************************************************
 * @file FaceData.cpp
 * @brief Implementation of face-centered data containers
 *****************************************************************************/

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "FaceData.hpp"

// ****************************** Constructor *******************************

template<typename T>
FaceData<T>::FaceData
(
    const std::string& fieldName,
    size_t numFaces
) : name_(fieldName),
    allFacesValues_(numFaces) {}

template<typename T>
FaceData<T>::FaceData
(
    const std::string& fieldName,
    size_t numFaces,
    const T& initialValue
) : name_(fieldName),
    allFacesValues_(numFaces, initialValue) {}

template<typename T>
T& FaceData<T>::operator[](size_t faceIndex)
{
    if (faceIndex >= allFacesValues_.size())
    {
        throw
            std::out_of_range
            (
                "Face index out of range in FaceField '"
              + name_ + "'"
            );
    }

    return allFacesValues_[faceIndex];
}

template<typename T>
const T& FaceData<T>::operator[](size_t faceIndex) const
{
    if (faceIndex >= allFacesValues_.size())
    {
        throw
            std::out_of_range
            (
                "Face index out of range in FaceField '"
              + name_ + "'"
            );
    }

    return allFacesValues_[faceIndex];
}

// ***************************** Setter Methods *****************************

template<typename T>
void FaceData<T>::setAll(const T& value)
{
    std::fill(allFacesValues_.begin(), allFacesValues_.end(), value);
}

// ***************************** Utility Methods ****************************

template<typename T>
void FaceData<T>::printSummary(size_t itemsToShow) const
{
    std::cout
        << "FaceField: " << name_ << " (Size: " << allFacesValues_.size()
        << ")" << std::endl;

    for (size_t i = 0; i < std::min(allFacesValues_.size(), itemsToShow); ++i)
    {
        std::cout
            << "  Face " << i << ": " << allFacesValues_[i] << std::endl;
    }

    if (allFacesValues_.size() > itemsToShow)
    {
        std::cout
            << "  ..." << std::endl;
    }
}

// ************************* Template Instantiations ***************************

template class FaceData<Scalar>;
template class FaceData<Vector>;
