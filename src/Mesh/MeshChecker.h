/******************************************************************************
 * @file MeshChecker.h
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

#include "Scalar.h"
#include "Mesh.h"


class MeshChecker
{
public:

    /// Construct a MeshChecker with a mesh view
    explicit MeshChecker(const Mesh& mesh) noexcept
    :
        mesh_(mesh)
    {}

    /// Copy constructor and assignment - Not copyable (const T& member)
    MeshChecker(const MeshChecker&) = delete;
    MeshChecker& operator=(const MeshChecker&) = delete;

    /// Move constructor and assignment - Not movable (const T& member)
    MeshChecker(MeshChecker&&) = delete;
    MeshChecker& operator=(MeshChecker&&) = delete;

    /// Destructor
    ~MeshChecker() noexcept = default;

    /// Perform mesh quality checks and report statistics
    void check() const;

private:

// Mesh quality thresholds

    /// Minimum face area before warning
    static constexpr Scalar minArea_ = S(1e-12);

    /// Minimum cell volume before warning
    static constexpr Scalar minVolume_ = S(1e-30);

    /// Non-orthogonality warning threshold (degrees)
    static constexpr Scalar maxNonOrthThreshold_ = S(70.0);

    /// Skewness warning threshold
    static constexpr Scalar maxSkewThreshold_ = S(4.0);

    /// Aspect ratio warning threshold
    static constexpr Scalar maxAspectThreshold_ = S(100.0);

// Private members

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

// Private methods

    /// Calculate face non-orthogonality using the OpenFOAM method
    [[nodiscard]] Scalar faceOrthogonality
    (
        const Vector& ownerCellCentroid,
        const Vector& neighborCellCentroid,
        const Vector& faceNormal
    ) const noexcept;

    /// Calculate face skewness using the OpenFOAM method
    [[nodiscard]] Scalar faceSkewness
    (
        const Face& face,
        const Vector& ownerCellCentroid,
        const Vector& neighborCellCentroid
    ) const;

    /// Calculate boundary face skewness
    [[nodiscard]] Scalar boundaryFaceSkewness
    (
        const Face& face,
        const Vector& ownerCellCentroid
    ) const;

    /// Calculate cell aspect ratio using the OpenFOAM method
    [[nodiscard]] Scalar cellAspectRatio(const Cell& cell) const;

    /// Validate mesh connectivity indices are in range
    bool validateConnectivity() const;

    /// Print up to 10 IDs from a list, with truncation
    static void printIndicesList
    (
        std::span<const size_t> indices,
        std::string_view entityName
    );
};
