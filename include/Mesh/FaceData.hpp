/******************************************************************************
 * @file FaceData.hpp
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

#include <algorithm>
#include <cstddef>
#include <span>
#include <vector>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Mesh.hpp"


template<typename T>
class FaceData
{
public:

    /// Construct zero-initialized field
    FaceData()
        : internalField_(Mesh::faceCount(), T{}) {}

    /**
     * @brief Construct field with initial value
     * @param initialValue Value to initialize all faces with
     */
    explicit FaceData(const T& initialValue)
        : internalField_(Mesh::faceCount(), initialValue) {}

// Setter methods

    /**
     * @brief Set all field values to a given value
     * @param value Value to assign to all faces
     */
    void setAll(const T& value)
    {
        std::fill(internalField_.begin(), internalField_.end(), value);
    }

// Accessor methods

    /**
     * @brief Get number of faces in the field
     * @return Number of faces
     */
    size_t size() const noexcept { return internalField_.size(); }

    /**
     * @brief Get pointer to field storage
     * @return Pointer to first element
     */
    T* data() noexcept { return internalField_.data(); }

    /**
     * @brief Get const pointer to field storage
     * @return Const pointer to first element
     */
    const T* data() const noexcept { return internalField_.data(); }

    /// Iterator access (range-based for loops)
    auto begin() noexcept { return internalField_.begin(); }
    auto end() noexcept { return internalField_.end(); }
    auto begin() const noexcept { return internalField_.begin(); }
    auto end() const noexcept { return internalField_.end(); }

// Operator methods

    /**
     * @brief Unchecked subscript operator
     * @param faceIndex Index of the face to access
     * @return Reference to field value at the face
     */
    T& operator[](size_t faceIndex) { return internalField_[faceIndex]; }

    /**
     * @brief Unchecked const subscript operator
     * @param faceIndex Index of the face to access
     * @return Const reference to field value at the face
     */
    const T& operator[](size_t faceIndex) const
    {
        return internalField_[faceIndex];
    }

// Helper methods

    /**
     * @brief Print field summary for debugging
     * @param itemsToShow Number of items to display
     */
    void printSummary(size_t itemsToShow) const;

private:

    /// Face-centered field values
    std::vector<T> internalField_;
};

/// Type alias for scalar face fields (e.g., mass flux)
using FaceFluxField = FaceData<Scalar>;

/// Type alias for vector face fields (e.g., gradients)
using FaceVectorField = FaceData<Vector>;
