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

#ifndef CELL_DATA_HPP
#define CELL_DATA_HPP

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

#include "Scalar.hpp"
#include "Vector.hpp"


template<typename T>
class CellData
{
public:

    /**
     * @brief Constructor for uninitialized field
     * @param fieldName Name identifier for the field
     * @param numCells Number of cells in the mesh
     */
    CellData
    (
        const std::string& fieldName,
        size_t numCells
    );

    /**
     * @brief Constructor with initial value
     * @param fieldName Name identifier for the field
     * @param numCells Number of cells in the mesh
     * @param initialValue Value to initialize all cells with
     */
    CellData
    (
        const std::string& fieldName,
        size_t numCells,
        const T& initialValue
    );

    /**
     * @brief Subscript operator
     * @param cellIndex Index of the cell to access
     * @return Reference to field value at the cell
     * @throws std::out_of_range if cellIndex is invalid
     */
    T& operator[](size_t cellIndex);

    /**
     * @brief Const subscript operator
     * @param cellIndex Index of the cell to access
     * @return Const reference to field value at the cell
     * @throws std::out_of_range if cellIndex is invalid
     */
    const T& operator[](size_t cellIndex) const;

    /**
     * @brief Get number of cells in the field
     * @return Number of cells
     */
    size_t size() const { return internalField_.size(); }

    /**
     * @brief Set all field values to a given value
     * @param value Value to assign to all cells
     */
    void setAll(const T& value);

    /**
     * @brief Get pointer to field storage
     * @return Pointer to first element
     */
    T* data() { return internalField_.data(); }

    /**
     * @brief Get const pointer to field storage
     * @return Const pointer to first element
     */
    const T* data() const { return internalField_.data(); }

    /**
     * @brief Get field name
     * @return Const reference to field name
     */
    const std::string& name() const { return name_; }

    /**
     * @brief Print field summary for debugging
     * @param itemsToShow Number of items to display
     */
    void printSummary(size_t itemsToShow) const;

private:

    /// Field name identifier
    std::string name_;

    /// Cell-centered field values
    std::vector<T> internalField_;
};

/// Type alias for velocity fields (Vector-valued)
using VelocityField = CellData<Vector>;

/// Type alias for pressure fields (Scalar-valued)
using PressureField = CellData<Scalar>;

/// Type alias for general scalar fields
using ScalarField = CellData<Scalar>;

/// Type alias for general vector fields
using VectorField = CellData<Vector>;

#endif // CELL_DATA_HPP
