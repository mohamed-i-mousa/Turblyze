/******************************************************************************
 * @file GradientScheme.hpp
 * @brief Gradient computation schemes for finite volume discretization
 *
 * @details This header defines the GradientScheme class, which provides 
 * robust gradient reconstruction methods for scalar and vector fields on
 * unstructured finite volume meshes. The implementation uses a weighted
 * least-squares approach with inverse distance squared weighting.
 *
 * @class GradientScheme
 *
 * The GradientScheme class provides:
 * - Cell-centered gradient computation via weighted least-squares
 * - Face-centered gradient interpolation with orthogonal correction
 * - Boundary condition handling for gradient stencils
 * - Distance-weighted averaging for internal face gradients
 *****************************************************************************/

#pragma once

#include <vector>
#include <array>
#include <string>
#include <optional>

#include "Scalar.hpp"
#include "Mesh.hpp"
#include "BoundaryConditions.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"


class GradientScheme
{
public:

    /**
     * @brief Construct gradient scheme with mesh context
     * @param mesh Mesh view (nodes, faces, cells)
     * @param bc Reference to boundary conditions manager
     */
    GradientScheme
    (
        const Mesh& mesh,
        const BoundaryConditions& bc
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    GradientScheme(const GradientScheme&) = delete;
    GradientScheme& operator=(const GradientScheme&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    GradientScheme(GradientScheme&&) = delete;
    GradientScheme& operator=(GradientScheme&&) = delete;

    /// Destructor
    ~GradientScheme() noexcept = default;

    /**
     * @brief Calculate gradient at a single cell using least-squares.
     * @param fieldName Name of the field for BC lookup
     * @param phi Scalar field for gradient calculation
     * @param cellIdx Index of the cell to compute gradient for
     * @param boundaryFaceValues Optional per-face boundary value overrides.
     *        The caller MUST initialise every entry to
     *        std::numeric_limits<Scalar>::quiet_NaN() before filling in
     *        specific overrides. Faces left as NaN fall back to the BC
     *        manager. Pass nullptr to bypass the override entirely
     * @param componentIdx For vector field components (0=x, 1=y, 2=z),
     *        used to extract scalar from vector BCs
     * @return Gradient vector at the specified cell
     */
    [[nodiscard]] Vector cellGradient
    (
        const std::string& fieldName,
        const ScalarField& phi,
        size_t cellIdx,
        const FaceData<Scalar>* boundaryFaceValues = nullptr,
        std::optional<int> componentIdx = std::nullopt
    ) const;

    /**
     * @brief Interpolate gradient at a single face
     * 
     * @details
     * For internal faces, computes a corrected face gradient using
     * distance-weighted interpolation of cell gradients with an
     * orthogonal correction to ensure consistency with the direct
     * cell-to-cell difference.
     * 
     * @param fieldName Name of the field for BC lookup
     * @param phi Cell-centered scalar field
     * @param gradPhiP Gradient at the owner cell
     * @param gradPhiN Gradient at the neighbor cell
     * @param faceIndex Index of the face
     * @return Gradient vector at the specified face
     */
    [[nodiscard]] Vector faceGradient
    (
        const std::string& fieldName,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        size_t faceIndex,
        std::optional<int> componentIdx = std::nullopt
    ) const;

    /**
     * @brief Apply Barth-Jespersen cell-based gradient limiter
     *
     * @details Scales each cell's gradient so that face-extrapolated
     * values do not exceed the range of neighboring cell values,
     * including values imposed by boundary conditions on boundary faces.
     * Prevents overshoot from steep gradients (e.g. near wall cells).
     *
     * @param fieldName Name of the field for BC lookup on boundary faces
     * @param phi Cell-centered scalar field
     * @param gradPhi Cell-centered gradient field (modified in place)
     * @param componentIdx For vector field components (0=x, 1=y, 2=z)
     */
    void limitGradient
    (
        const std::string& fieldName,
        const ScalarField& phi,
        VectorField& gradPhi,
        std::optional<int> componentIdx = std::nullopt
    ) const;


private:

// Private methods

    /**
     * @brief Distance-weighted linear interpolation of gradients
     * @param face Internal face for interpolation
     * @param gradPhiP Gradient at the owner cell
     * @param gradPhiN Gradient at the neighbor cell
     * @return Linearly interpolated gradient vector at the face
     */
    Vector averageFaceGradient
    (
        const Face& face,
        const Vector& gradPhiP,
        const Vector& gradPhiN
    ) const;

    /**
     * @brief Calculate boundary face gradient based on BC type
     * @param face Boundary face
     * @param fieldName Name of the field for BC lookup
     * @param phi Scalar field values
     * @param cellGradient Gradient at the owner cell
     * @return Gradient vector at the boundary face
     */
    Vector boundaryFaceGradient
    (
        const std::string& fieldName,
        const ScalarField& phi,
        const Vector& cellGradient,
        const Face& face,
        std::optional<int> componentIdx = std::nullopt
    ) const;

    /// Pre-compute inverse of ATA matrix for each cell
    void precomputeInverseATA();

// Private members

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

    /// Reference to boundary conditions manager
    const BoundaryConditions& bcManager_;

    /// Cached inverse of ATA per cell {xx, xy, xz, yy, yz, zz}
    std::vector<std::array<Scalar, 6>> invATA_;

    /// Minimum fraction of ||dPf|| used as normal distance to a boundary face.
    /// Prevents gradient amplification beyond ~87 degrees of non-orthogonality.
    static constexpr Scalar minNormalFraction_ = S(0.05);
};
