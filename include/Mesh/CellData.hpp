/******************************************************************************
 * @file CellData.hpp
 * @brief Template container for cell-centered field data storage
 * 
 * @details This header defines a generic template class for storing field 
 * variables at cell centers in finite volume meshes. The container manages 
 * cell-based data including velocity fields, pressure field, and 
 * cell-centered gradients
 * 
 * @class CellData<T>
 * @tparam T Type of field value stored at each cell (e.g., Scalar, Vector)
 *
 * The CellData template provides:
 * - Type-safe storage for cell-centered field variables
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
#include "Tensor.hpp"
#include "Mesh.hpp"


template<typename T>
class CellData
{
public:

    /// Construct zero-initialized field
    CellData()
    :
        internalField_(Mesh::cellCount(), T{})
    {}

    /**
     * @brief Construct field with initial value
     * @param initialValue Value to initialize all cells with
     */
    explicit CellData(const T& initialValue)
    :
        internalField_(Mesh::cellCount(), initialValue)
    {}
    
// Setter methods

    /**
     * @brief Set all field values to a given value
     * @param value Value to assign to all cells
     */
    void setAll(const T& value)
    {
        std::fill(internalField_.begin(), internalField_.end(), value);
    }

// Accessor methods

    /**
     * @brief Get number of cells in the field
     * @return Number of cells
     */
    [[nodiscard]] size_t size() const noexcept
    {
        return internalField_.size();
    }

    /**
     * @brief Get pointer to field storage
     * @return Pointer to first element
     */
    [[nodiscard]] T* data() noexcept { return internalField_.data(); }

    /**
     * @brief Get const pointer to field storage
     * @return Const pointer to first element
     */
    [[nodiscard]] const T* data() const noexcept
    {
        return internalField_.data();
    }

    /// Iterator access (range-based for loops)
    [[nodiscard]] auto begin() noexcept { return internalField_.begin(); }
    [[nodiscard]] auto end() noexcept { return internalField_.end(); }
    [[nodiscard]] auto begin() const noexcept { return internalField_.begin(); }
    [[nodiscard]] auto end() const noexcept { return internalField_.end(); }

// Operator methods

    /**
     * @brief Unchecked subscript operator
     * @param cellIndex Index of the cell to access
     * @return Reference to field value at the cell
     */
    T& operator[](size_t cellIndex) { return internalField_[cellIndex]; }

    /**
     * @brief Unchecked const subscript operator
     * @param cellIndex Index of the cell to access
     * @return Const reference to field value at the cell
     */
    const T& operator[](size_t cellIndex) const
    {
        return internalField_[cellIndex];
    }

// Helper methods

    /**
     * @brief Print field summary for debugging
     * @param itemsToShow Number of items to display
     */
    void printSummary(size_t itemsToShow) const;

private:

    /// Cell-centered field values
    std::vector<T> internalField_;
};

/// Type alias for general scalar fields
using ScalarField = CellData<Scalar>;

/// Type alias for general vector fields
using VectorField = CellData<Vector>;

/// Type alias for general tensor fields
using TensorField = CellData<Tensor>;
