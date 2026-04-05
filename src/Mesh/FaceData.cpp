/******************************************************************************
 * @file FaceData.cpp
 * @brief Implementation of face-centered data containers
 *****************************************************************************/

#include "FaceData.hpp"

#include <algorithm>
#include <iostream>
#include <utility>

#include "ErrorHandler.hpp"

// ****************************** Constructor *******************************

template<typename T>
FaceData<T>::FaceData
(
    std::string fieldName,
    size_t numFaces
) : name_(std::move(fieldName)),
    internalField_(numFaces)
{
    if (numFaces == 0)
    {
        FatalError("FaceData: cannot create field with zero faces");
    }
}

template<typename T>
FaceData<T>::FaceData
(
    std::string fieldName,
    size_t numFaces,
    const T& initialValue
) : name_(std::move(fieldName)),
    internalField_(numFaces, initialValue)
{
    if (numFaces == 0)
    {
        FatalError("FaceData: cannot create field with zero faces");
    }
}

// ***************************** Setter Methods *****************************

template<typename T>
void FaceData<T>::setAll(const T& value)
{
    std::fill(internalField_.begin(), internalField_.end(), value);
}

// ***************************** Utility Methods ****************************

template<typename T>
void FaceData<T>::printSummary(size_t itemsToShow) const
{
    std::cout
        << "FaceData: " << name_ << " (Size: " << internalField_.size()
        << ")\n";

    for (size_t i = 0; i < std::min(internalField_.size(), itemsToShow); ++i)
    {
        std::cout
            << "  Face " << i << ": " << internalField_[i] << '\n';
    }

    if (internalField_.size() > itemsToShow)
    {
        std::cout
            << "  ...\n";
    }
}

// ************************* Template Instantiations ***************************

template class FaceData<Scalar>;
template class FaceData<Vector>;
