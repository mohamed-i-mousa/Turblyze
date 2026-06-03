/******************************************************************************
 * @file FaceData.h
 * @brief Template container for face-centered field data storage
 *
 * @details This header defines a generic template class for storing field
 * variables at face centers in finite volume meshes. The container manages
 * face-based data including mass fluxes, face velocities, and face-centered
 * gradients
 *
 * @class FaceData<T>
 * @tparam T Type of field value stored at each face (e.g., Scalar, Vector)
 *
 * The FaceData template provides:
 * - Type-safe storage for face-centered field variables
 * - Face-specific initialization and assignment operations
 * - Debugging output for face field analysis
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <algorithm>
#include <span>
#include <vector>
#include <concepts>
#include <iostream>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "Mesh.h"
#include "Integer.h"

// *************************** concept FaceFieldType **************************

template<typename T>
concept FaceFieldType = std::same_as<T, Scalar> || std::same_as<T, Vector>;

// ****************************** class FaceData ******************************

template<FaceFieldType T>
class FaceData
{
public:

// ************************* Special Member Functions *************************

    /// Construct zero-initialized field
    FaceData()
    :
        internalField_(Mesh::faceCount(), T{})
    {}

    /// Construct field with initial value
    explicit FaceData(const T& initialValue)
    :
        internalField_(Mesh::faceCount(), initialValue)
    {}

// ****************************** Setter Methods ******************************

    /// Set all field values to a given value
    void setAll(const T& value)
    {
        std::fill(internalField_.begin(), internalField_.end(), value);
    }

// ***************************** Accessor Methods *****************************

    /// Get number of faces in the field
    [[nodiscard]] Count size() const noexcept
    {
        return internalField_.size();
    }

    /// Get pointer to field storage
    [[nodiscard]] T* data() noexcept
    {
        return internalField_.data();
    }

    /// Get const pointer to field storage
    [[nodiscard]] const T* data() const noexcept
    {
        return internalField_.data();
    }

    /// Iterator access (range-based for loops)
    [[nodiscard]] auto begin() noexcept
    {
        return internalField_.begin();
    }

    [[nodiscard]] auto end() noexcept
    {
        return internalField_.end();
    }

    [[nodiscard]] auto begin() const noexcept
    {
        return internalField_.begin();
    }

    [[nodiscard]] auto end() const noexcept
    {
        return internalField_.end();
    }

// ***************************** Operator Methods *****************************

    /// Unchecked subscript operator
    T& operator[](Index faceIndex) noexcept
    {
        return internalField_[faceIndex];
    }

    /// Unchecked const subscript operator
    const T& operator[](Index faceIndex) const noexcept
    {
        return internalField_[faceIndex];
    }

// ****************************** Helper Methods ******************************

    /// Print field summary for debugging
    void printSummary(Count itemsToShow) const
    {
        std::cout
            << "FaceData (Size: " << internalField_.size() << ")\n";

        const Count count = std::min(internalField_.size(), itemsToShow);

        for
        (
            Index faceIdx = 0;
            faceIdx < count;
            ++faceIdx
        )
        {
            std::cout
                << "  Face " << faceIdx << ": " << internalField_[faceIdx]
                << '\n';
        }

        if (internalField_.size() > itemsToShow)
        {
            std::cout
                << "  ...\n";
        }
    }

// ****************************** Private Members *****************************

private:

    /// Face-centered field values
    std::vector<T> internalField_;
};

// ********************************** Aliases *********************************

/// Type alias for scalar face fields (e.g., mass flux)
using FaceFluxField = FaceData<Scalar>;

/// Type alias for vector face fields (e.g., gradients)
using FaceVectorField = FaceData<Vector>;
