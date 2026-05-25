/******************************************************************************
 * @file GradientScheme.h
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

#include "Scalar.h"
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "Field.h"


class GradientScheme
{
public:

    /// Construct gradient scheme with mesh context
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

    /// Calculate gradient at a single cell using least-squares
    [[nodiscard]] Vector cellGradient
    (
        Field field,
        const ScalarField& phi,
        size_t cellIdx
    ) const;

    /// Interpolate gradient at a single face
    [[nodiscard]] Vector faceGradient
    (
        Field field,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        size_t faceIndex
    ) const;

    /// Apply cell-based gradient limiter
    void limitGradient
    (
        Field field,
        const ScalarField& phi,
        VectorField& gradPhi
    ) const;

    /// Compute a limited cell-centered gradient field
    void fieldGradient
    (
        Field field,
        const ScalarField& phi,
        VectorField& gradPhi
    ) const;


private:

// Private methods

    /// Distance-weighted linear interpolation of gradients
    Vector averageFaceGradient
    (
        const Face& face,
        const Vector& gradPhiP,
        const Vector& gradPhiN
    ) const;

    /// Calculate boundary face gradient based on BC type
    Vector boundaryFaceGradient
    (
        Field field,
        const ScalarField& phi,
        const Vector& cellGradient,
        const Face& face
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
