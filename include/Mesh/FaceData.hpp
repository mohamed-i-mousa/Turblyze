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

#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include "Scalar.hpp"
#include "Vector.hpp"


template<typename T>
class FaceData
{
public:

    /// Default constructor
    FaceData() = default;

    /**
     * @brief Constructor for uninitialized field
     * @param fieldName Name identifier for the field
     * @param numFaces Number of faces in the mesh
     */
    FaceData
    (
        std::string fieldName,
        size_t numFaces
    );

    /**
     * @brief Constructor with initial value
     * @param fieldName Name identifier for the field
     * @param numFaces Number of faces in the mesh
     * @param initialValue Value to initialize all faces with
     */
    FaceData
    (
        std::string fieldName,
        size_t numFaces,
        const T& initialValue
    );

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

    /**
     * @brief Get number of faces in the field
     * @return Number of faces
     */
    size_t size() const noexcept { return internalField_.size(); }

    /**
     * @brief Set all field values to a given value
     * @param value Value to assign to all faces
     */
    void setAll(const T& value);

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

    /**
     * @brief Get a mutable view of the field storage
     * @return Non-owning span over face values
     */
    std::span<T> span() noexcept { return internalField_; }

    /**
     * @brief Get a read-only view of the field storage
     * @return Non-owning span over const face values
     */
    std::span<const T> span() const noexcept { return internalField_; }

    /// Iterator access (range-based for loops)
    auto begin() noexcept { return internalField_.begin(); }
    auto end() noexcept { return internalField_.end(); }
    auto begin() const noexcept { return internalField_.begin(); }
    auto end() const noexcept { return internalField_.end(); }

    /**
     * @brief Get field name
     * @return Const reference to field name
     */
    const std::string& name() const noexcept { return name_; }

    /**
     * @brief Print field summary for debugging
     * @param itemsToShow Number of items to display
     */
    void printSummary(size_t itemsToShow) const;

private:

    /// Field name identifier
    std::string name_;

    /// Face-centered field values
    std::vector<T> internalField_;
};

/// Type alias for scalar face fields (e.g., mass flux)
using FaceFluxField = FaceData<Scalar>;

/// Type alias for vector face fields (e.g., gradients)
using FaceVectorField = FaceData<Vector>;
