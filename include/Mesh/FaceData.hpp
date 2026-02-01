/******************************************************************************
 * @file FaceData.hpp
 * @brief Template container for face-centered field data storage
 * 
 * This header defines a generic template class for storing field variables
 * at face centers in finite volume meshes. The container manages face-based
 * data including mass fluxes, face velocities, and face-centered gradients
 * 
 * @class FaceData<T>
 * @tparam T Type of field value stored at each face (e.g., Scalar, Vector)
 *
 * The FaceData template provides:
 * - Type-safe storage for face-centered field variables
 * - Face-specific initialization and assignment operations  
 * - Debugging output for face field analysis
 *****************************************************************************/

#ifndef FACE_DATA_HPP
#define FACE_DATA_HPP

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

#include "Scalar.hpp"
#include "Vector.hpp"


template<typename T>
class FaceData
{
public:

    /**
     * @brief Constructor for uninitialized field
     * @param fieldName Name identifier for the field
     * @param numFaces Number of faces in the mesh
     */
    FaceData
    (
        const std::string& fieldName,
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
        const std::string& fieldName,
        size_t numFaces,
        const T& initialValue
    );

    /**
     * @brief Subscript operator
     * @param faceIndex Index of the face to access
     * @return Reference to field value at the face
     * @throws std::out_of_range if faceIndex is invalid
     */
    T& operator[](size_t faceIndex);

    /**
     * @brief Const subscript operator
     * @param faceIndex Index of the face to access
     * @return Const reference to field value at the face
     * @throws std::out_of_range if faceIndex is invalid
     */
    const T& operator[](size_t faceIndex) const;

    /**
     * @brief Get number of faces in the field
     * @return Number of faces
     */
    size_t size() const { return allFacesValues_.size(); }

    /**
     * @brief Set all field values to a given value
     * @param value Value to assign to all faces
     */
    void setAll(const T& value);

    /**
     * @brief Print field summary for debugging
     * @param itemsToShow Number of items to display
     */
    void printSummary(size_t itemsToShow) const;

    /**
     * @brief Get field name
     * @return Const reference to field name
     */
    const std::string& getName() const { return name_; }

    /**
     * @brief Set field name
     * @param fieldName New name for the field
     */
    void setName(const std::string& fieldName) { name_ = fieldName; }

    /**
     * @brief Get all faces values (const)
     * @return Const reference to face values vector
     */
    const std::vector<T>& getAllFacesValues() const { return allFacesValues_; }

    /**
     * @brief Get all faces values (mutable)
     * @return Reference to face values vector
     */
    std::vector<T>& getAllFacesValues() { return allFacesValues_; }

private:

    /// Field name identifier
    std::string name_;

    /// Face-centered field values
    std::vector<T> allFacesValues_;
};

/// Type alias for scalar face fields (e.g., mass flux)
using FaceFluxField = FaceData<Scalar>;

/// Type alias for vector face fields (e.g., gradients)
using FaceVectorField = FaceData<Vector>;

#endif // FACE_DATA_HPP