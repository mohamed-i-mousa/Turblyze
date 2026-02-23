/******************************************************************************
 * @file MeshChecker.hpp
 * @brief Mesh quality assessment and diagnostic utilities
 *
 * @details This header provides mesh quality checking functions that analyze 
 * geometric properties and report statistics to help identify potential 
 * numerical issues.
 *
 * @class MeshChecker
 *
 * Key mesh quality metrics:
 * - Minimum and maximum face areas
 * - Minimum and maximum cell volumes
 * - Skewness
 * - Non-Orthogonality
 * - Aspect Ratio
 *****************************************************************************/

#ifndef MESH_CHECKER_HPP
#define MESH_CHECKER_HPP

#include <vector>
#include "Face.hpp"
#include "Cell.hpp"
#include "Vector.hpp"
#include "Scalar.hpp"


class MeshChecker
{
public:

    /**
     * @brief Construct a MeshChecker with mesh data references
     * @param nodes Vector containing all mesh nodes
     * @param faces Vector containing all mesh faces
     * @param cells Vector containing all mesh cells
     */
    MeshChecker
    (
        const std::vector<Vector>& nodes,
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells
    );

    /// Perform mesh quality checks and report statistics
    void check() const;

private:

// Mesh quality thresholds

    /// Minimum face area before warning
    static constexpr Scalar minArea_ = 1e-12;

    /// Minimum cell volume before warning
    static constexpr Scalar minVolume_ = 1e-30;

// Private members

    /// Reference to all mesh nodes
    const std::vector<Vector>& allNodes_;

    /// Reference to all mesh faces
    const std::vector<Face>& allFaces_;

    /// Reference to all mesh cells
    const std::vector<Cell>& allCells_;

// Private methods

    /**
     * @brief Calculate face non-orthogonality (OpenFOAM method)
     * @param ownerCellCentroid Owner cell center
     * @param neighborCellCentroid Neighbor cell center
     * @param faceNormal Face normal vector (unit vector)
     * @return Cosine of angle between cell centers line and face normal
     */
    Scalar calculateFaceOrthogonality
    (
        const Vector& ownerCellCentroid,
        const Vector& neighborCellCentroid,
        const Vector& faceNormal
    ) const;

    /**
     * @brief Calculate face skewness (OpenFOAM method)
     * @param face Face object for accessing vertices
     * @param ownerCellCentroid Owner cell center
     * @param neighborCellCentroid Neighbor cell center
     * @param faceCentroid Face center
     * @param faceNormal Face normal vector
     * @return Skewness value
     */
    Scalar calculateFaceSkewness
    (
        const Face& face,
        const Vector& ownerCellCentroid,
        const Vector& neighborCellCentroid,
        const Vector& faceCentroid,
        const Vector& faceNormal
    ) const;

    /**
     * @brief Calculate boundary face skewness
     * @param face Face object for accessing vertices
     * @param ownerCellCentroid Owner cell center
     * @param faceCentroid Face center
     * @param faceNormal Face normal vector
     * @return Boundary skewness value
     */
    Scalar calculateBoundarySkewness
    (
        const Face& face,
        const Vector& ownerCellCentroid,
        const Vector& faceCentroid,
        const Vector& faceNormal
    ) const;

    /**
     * @brief Calculate cell aspect ratio (OpenFOAM method)
     * @param cell Cell object
     * @return Aspect ratio (1.0 = perfect cube, higher = elongated)
     */
    Scalar calculateCellAspectRatio(const Cell& cell) const;
};

#endif // MESH_CHECKER_HPP
