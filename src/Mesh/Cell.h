/******************************************************************************
 * @file Cell.h
 * @brief Represents a computational cell in the mesh
 *
 * @details This header defines the Cell class, which represents a finite
 * control volume in the computational mesh. The cell is the primary entity
 * where flow variables (pressure, velocity, etc.) are stored and solved.
 * The cell is defined by a collection of bounding faces that form a closed
 * volume.
 *
 * @class Cell
 * - Topological connectivity (bounding faces, neighboring cells)
 * - Orientation data (face alignment signs relative to cell)
 * - Cell properties (centroid, volume)
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <cstdint>
#include <vector>
#include <span>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "Face.h"

// ******************************** class Cell ********************************

class Cell
{
public:

// ************************* Special Member Functions *************************

    /// Default constructor
    Cell() noexcept = default;

    /**
     * @brief Constructs cell with connectivity data
     * @param cellIdx Unique cell identifier
     * @param faces Indices of bounding faces
     * @param neighbors Indices of neighboring cells
     * @param signs Face normal direction signs
     */
    Cell
    (
        size_t cellIdx,
        std::vector<size_t> faces,
        std::vector<size_t> neighbors,
        std::vector<int8_t> signs
    )
    : 
        idx_(cellIdx),
        faceIndices_(std::move(faces)),
        neighborCellIndices_(std::move(neighbors)),
        faceSigns_(std::move(signs))
    {}

// ****************************** Setter Methods ******************************

    /// Set cell identifier
    void setIdx(size_t cellIdx) noexcept
    {
        idx_ = cellIdx;
    }

    /// Add a bounding face with its normal direction sign
    void addFace(size_t faceIdx, int8_t sign)
    {
        faceIndices_.push_back(faceIdx);
        faceSigns_.push_back(sign);
    }

    /// Set all neighbor cell indices
    void setNeighborCellIndices(std::span<const size_t> neighbors)
    {
        neighborCellIndices_.assign(neighbors.begin(), neighbors.end());
    }

// ***************************** Accessor Methods *****************************

    /// Get cell identifier
    [[nodiscard]] size_t idx() const noexcept
    {
        return idx_;
    }

    /// Get bounding face indices
    [[nodiscard]] std::span<const size_t> faceIndices() const noexcept
    {
        return faceIndices_;
    }

    /// Get neighboring cell indices
    [[nodiscard]] std::span<const size_t> neighborCellIndices() const noexcept
    {
        return neighborCellIndices_;
    }

    /// Get face normal direction signs
    [[nodiscard]] std::span<const int8_t> faceSigns() const noexcept
    {
        return faceSigns_;
    }

    /// Get cell centroid
    [[nodiscard]] const Vector& centroid() const noexcept
    {
        return centroid_;
    }

    /// Get cell volume
    [[nodiscard]] Scalar volume() const noexcept
    {
        return volume_;
    }

// ************************ Geometric Property Methods ************************

    /// Check if geometric properties calculated
    [[nodiscard]] bool geometricPropertiesCalculated() const noexcept
    {
        return geometricPropertiesCalculated_;
    }

    /// Calculate cell volume and centroid
    void geometricProperties
    (
        std::span<const Face> allFaces,
        std::span<const FaceIntegrals> allFaceIntegrals
    );

// ****************************** Private Members *****************************

private:

    /// Unique cell identifier
    size_t idx_ = 0;

    /// Indices of faces that bound this cell
    std::vector<size_t> faceIndices_;

    /// Indices of neighboring cells
    std::vector<size_t> neighborCellIndices_;

    /// Face normal direction signs (+1 outward, -1 inward)
    std::vector<int8_t> faceSigns_;

    /// Cell geometric center
    Vector centroid_;

    /// Cell volume
    Scalar volume_ = S(0.0);

    /// Flag indicating if geometry has been calculated
    bool geometricPropertiesCalculated_ = false;

};

// *************************** Non-Member Functions ***************************

/// Stream output operator for Cell
std::ostream& operator<<(std::ostream& os, const Cell& c);
