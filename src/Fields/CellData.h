/******************************************************************************
 * @file CellData.h
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
#include "Tensor.h"
#include "Mesh.h"
#include "Integer.h"

// ************************* concept CellFieldType ****************************

template<typename T>
concept CellFieldType =
    std::same_as<T, Scalar>
 || std::same_as<T, Vector>
 || std::same_as<T, Tensor>;

// ****************************** class CellData ******************************

template<CellFieldType T>
class CellData
{
public:

// ************************* Special Member Functions *************************

    /// Construct zero-initialized field
    CellData()
    :
        internalField_(Mesh::cellCount(), T{})
    {}

    /// Construct field with initial value
    explicit CellData(const T& initialValue)
    :
        internalField_(Mesh::cellCount(), initialValue)
    {}

// ****************************** Setter Methods ******************************

    /// Set all field values to a given value
    void setAll(const T& value)
    {
        std::fill(internalField_.begin(), internalField_.end(), value);
    }

// ***************************** Accessor Methods *****************************

    /// Get number of cells in the field
    [[nodiscard]] Count size() const noexcept
    {
        return internalField_.size();
    }

    /// Check whether the field contains no cells
    [[nodiscard]] bool empty() const noexcept
    {
        return internalField_.empty();
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
    T& operator[](Index cellIndex) noexcept
    {
        return internalField_[cellIndex];
    }

    /// Unchecked const subscript operator
    const T& operator[](Index cellIndex) const noexcept
    {
        return internalField_[cellIndex];
    }

// ****************************** Helper Methods ******************************

    /// Print field summary for debugging
    void printSummary(Count itemsToShow) const
    {
        std::cout
            << "CellData (Size: " << internalField_.size() << ")\n";

        const Count count = std::min(internalField_.size(), itemsToShow);

        for
        (
            Index cellIdx = 0;
            cellIdx < count;
            ++cellIdx
        )
        {
            std::cout
                << "  Cell " << cellIdx << ": " << internalField_[cellIdx]
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

    /// Cell-centered field values
    std::vector<T> internalField_;
};

// ********************************** Aliases *********************************

/// Type alias for general scalar fields
using ScalarField = CellData<Scalar>;

/// Type alias for general vector fields
using VectorField = CellData<Vector>;

/// Type alias for general tensor fields
using TensorField = CellData<Tensor>;
