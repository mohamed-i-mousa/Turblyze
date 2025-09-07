/******************************************************************************
 * @file FaceData.h
 * @brief Template container for face-centered field data storage
 * 
 * This header defines a generic template class for storing field variables
 * at face centers in finite volume meshes. The container manages face-based
 * data including mass fluxes, face velocities, and face-centered gradients
 * with type safety and bounds checking. It handles the specialized storage
 * requirements for face-based quantities in collocated finite volume schemes.
 * 
 * @class FaceData<T>
 * 
 * The FaceData template provides:
 * - Type-safe storage for face-centered field variables
 * - Bounds-checked access with clear error reporting
 * - Face-specific initialization and assignment operations  
 * - Integration with Rhie-Chow interpolation procedures
 * - Support for mass flux conservation and face velocity storage
 * - Debugging output capabilities for face field analysis
 * 
 * Common instantiations:
 * - FaceFluxField = FaceData<Scalar> for mass fluxes, flow rates
 * - FaceVectorField = FaceData<Vector> for face velocities, face gradients
 *****************************************************************************/

#ifndef FACEDATA_H
#define FACEDATA_H

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

#include "Scalar.h"
#include "Vector.h"

/**
 * @brief Template class for storing face-centered field data
 * @tparam T Data type (Scalar, Vector)
 * 
 * Provides storage and access methods for field data defined at face centers.
 * Commonly used for fluxes, gradients, and other face-based quantities
 * in finite volume calculations.
 */
template<typename T>
class FaceData 
{
public:

    /// Field name identifier
    std::string name;
    
    /// Face-centered field values
    std::vector<T> allFacesValues;

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
    inline size_t size() const { return allFacesValues.size(); }
    
    /**
     * @brief Set all field values to a given value
     * @param value Value to assign to all faces
     */
    void setAll(const T& value);
    
    /**
     * @brief Print field summary for debugging
     * @param itemsToShow Number of items to display (default: 5)
     */
    void printSummary(size_t itemsToShow = 5) const;
};

/// Type alias for scalar face fields (e.g., mass flux)
using FaceFluxField = FaceData<Scalar>;

/// Type alias for vector face fields (e.g., gradients)
using FaceVectorField = FaceData<Vector>;

#endif