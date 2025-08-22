#ifndef CELLDATA_H
#define CELLDATA_H

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

#include "Scalar.h"
#include "Vector.h"

/**
 * @brief Template class for storing cell-centered field data
 * @tparam T Data type (Scalar, Vector, etc.)
 * 
 * Provides storage and access methods for field data defined at cell centers.
 * Supports both scalar and vector fields with bounds checking and 
 * initialization options.
 */
template<typename T>
class CellData 
{
public:
    /// Field name identifier
    std::string name;
    
    /// Cell-centered field values
    std::vector<T> internalField;

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
    size_t size() const;
    
    /**
     * @brief Set all field values to a given value
     * @param value Value to assign to all cells
     */
    void setAll(const T& value);
    
    /**
     * @brief Print field summary for debugging
     * @param itemsToShow Number of items to display
     */
    void printSummary(size_t itemsToShow) const;
};

/// Type alias for velocity fields (Vector-valued)
using VelocityField = CellData<Vector>;

/// Type alias for pressure fields (Scalar-valued)
using PressureField = CellData<Scalar>;

/// Type alias for general scalar fields
using ScalarField = CellData<Scalar>;

/// Type alias for general vector fields
using VectorField = CellData<Vector>;

#endif