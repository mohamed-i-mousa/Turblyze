/******************************************************************************
 * @file FaceData.cpp
 * @brief Implementation of face-centered data containers
 *****************************************************************************/

#include "FaceData.h"

#include <iostream>

#include "ErrorHandler.h"

// ***************************** Helper Methods *****************************

template<typename T>
void FaceData<T>::printSummary(size_t itemsToShow) const
{
    std::cout
        << "FaceData (Size: " << internalField_.size() << ")\n";

    const size_t count = std::min(internalField_.size(), itemsToShow);

    for
    (
        size_t faceIdx = 0;
        faceIdx < count;
        ++faceIdx
    )
    {
        std::cout
            << "  Face " << faceIdx << ": " << internalField_[faceIdx] << '\n';
    }

    if (internalField_.size() > itemsToShow)
    {
        std::cout
            << "  ...\n";
    }
}

// ************************* Template Instantiations **************************

template class FaceData<Scalar>;
template class FaceData<Vector>;
