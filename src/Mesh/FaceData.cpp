#include "FaceData.h"
#include <stdexcept>
#include <iostream>
#include <algorithm>

template<typename T>
FaceData<T>::FaceData
(
    const std::string& fieldName,
    size_t numFaces
) : name(fieldName), 
    allFacesValues(numFaces) {}

template<typename T>
FaceData<T>::FaceData
(
    const std::string& fieldName,
    size_t numFaces,
    const T& initialValue
) : name(fieldName),
    allFacesValues(numFaces, initialValue) {}

template<typename T>
T& FaceData<T>::operator[](size_t faceIndex)
{
    if (faceIndex >= allFacesValues.size())
    {
        throw std::out_of_range
        (
            "Face index out of range in FaceField '" + name + "'"
        );
    }

    return allFacesValues[faceIndex];
}

template<typename T>
const T& FaceData<T>::operator[](size_t faceIndex) const
{
    if (faceIndex >= allFacesValues.size())
    {
        throw std::out_of_range
        (
            "Face index out of range in FaceField '" + name + "'"
        );
    }

    return allFacesValues[faceIndex];
}


template<typename T>
void FaceData<T>::setAll(const T& value)
{
    for (size_t i = 0; i < allFacesValues.size(); i++)
    {
        allFacesValues[i] = value;
    }
}

template<typename T>
void FaceData<T>::printSummary(size_t itemsToShow) const
{
    std::cout   << "FaceField: " << name << " (Size: " 
                << allFacesValues.size() << ")" << std::endl;

    for (size_t i = 0; i < std::min(allFacesValues.size(), itemsToShow); ++i)
    {
        std::cout   << "  Face " << i << ": " << allFacesValues[i] 
                    << std::endl;
    }
    
    if (allFacesValues.size() > itemsToShow)
    {
        std::cout   << "  ..." << std::endl;
    }
}

template class FaceData<Scalar>;
template class FaceData<Vector>;
