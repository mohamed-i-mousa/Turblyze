/******************************************************************************
 * @file MeshChecker.hpp
 * @brief Mesh quality assessment and diagnostic utilities
 *
 * @details This header provides mesh quality checking functions that analyze
 * geometric properties and report statistics to help identify potential
 * numerical issues.
 *
 * @class MeshChecker
 * Key mesh quality metrics:
 * - Minimum and maximum face areas
 * - Minimum and maximum cell volumes
 * - Skewness
 * - Non-Orthogonality
 * - Aspect Ratio
 *****************************************************************************/

#pragma once

#include <vector>

#include "Scalar.hpp"
#include "Mesh.hpp"


class MeshChecker
{
public:

    /**
     * @brief Construct a MeshChecker with a mesh view
     * @param mesh Mesh view (nodes, faces, cells)
     */
    MeshChecker(const Mesh& mesh) noexcept
    :
        mesh_(mesh)
    {}

    /// Perform mesh quality checks and report statistics
    void check() const;

private:

// Mesh quality thresholds

    /// Minimum face area before warning
    static constexpr Scalar minArea_ = 1e-12;

    /// Minimum cell volume before warning
    static constexpr Scalar minVolume_ = 1e-30;

    /// Non-orthogonality warning threshold (degrees)
    static constexpr Scalar maxNonOrthThreshold_ = 70.0;

    /// Skewness warning threshold
    static constexpr Scalar maxSkewThreshold_ = 4.0;

    /// Aspect ratio warning threshold
    static constexpr Scalar maxAspectThreshold_ = 100.0;

// Private members

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

// Private methods

    /**
     * @brief Calculate face non-orthogonality (OpenFOAM method)
     * @param ownerCellCentroid Owner cell center
     * @param neighborCellCentroid Neighbor cell center
     * @param faceNormal Face normal vector (unit vector)
     * @return Cosine of angle between cell centers line and face normal
     */
    [[nodiscard]] Scalar faceOrthogonality
    (
        const Vector& ownerCellCentroid,
        const Vector& neighborCellCentroid,
        const Vector& faceNormal
    ) const noexcept;

    /**
     * @brief Calculate face skewness (OpenFOAM method)
     * @param face Face object for accessing vertices and geometry
     * @param ownerCellCentroid Owner cell center
     * @param neighborCellCentroid Neighbor cell center
     * @return Skewness value
     */
    [[nodiscard]] Scalar faceSkewness
    (
        const Face& face,
        const Vector& ownerCellCentroid,
        const Vector& neighborCellCentroid
    ) const;

    /**
     * @brief Calculate boundary face skewness
     * @param face Face object for accessing vertices and geometry
     * @param ownerCellCentroid Owner cell center
     * @return Boundary skewness value
     */
    [[nodiscard]] Scalar boundaryFaceSkewness
    (
        const Face& face,
        const Vector& ownerCellCentroid
    ) const;

    /**
     * @brief Calculate cell aspect ratio (OpenFOAM method)
     * @param cell Cell object
     * @return Aspect ratio (1.0 = perfect cube, higher = elongated)
     */
    [[nodiscard]] Scalar cellAspectRatio(const Cell& cell) const;

    /**
     * @brief Validate mesh connectivity indices are in range
     * @return True if all indices are valid, false otherwise
     */
    bool validateConnectivity() const;

    /// Print up to 10 IDs from a list, with truncation
    static void printIndicesList
    (
        std::span<const size_t> indices,
        std::string_view entityName
    );
};
