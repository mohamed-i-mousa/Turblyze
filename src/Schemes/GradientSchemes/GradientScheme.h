/******************************************************************************
 * @file GradientScheme.h
 * @brief Abstract base class for finite volume gradient reconstruction
 *
 * @details This header defines the GradientScheme abstract base class, which
 * provides scheme-independent gradient services for scalar and vector fields
 * on unstructured finite volume meshes. The single scheme-specific operation,
 * the cell-centered gradient, is declared pure virtual so concrete schemes
 * (e.g. LeastSquares) can override it.
 *
 * @class GradientScheme
 * - Cell-centered gradient computation (pure virtual, scheme-specific)
 * - Face-centered gradient interpolation with orthogonal correction
 * - Cell-based gradient limiting (Barth-Jespersen)
 * - Distance-weighted averaging for internal face gradients
 * - Boundary condition handling for face gradients
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Scalar.h"
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "Field.h"

// *************************** class GradientScheme ***************************

class GradientScheme
{
public:

// ************************* Special Member Functions *************************

    /// Copy constructor and assignment - Not copyable (const T& members)
    GradientScheme(const GradientScheme&) = delete;
    GradientScheme& operator=(const GradientScheme&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    GradientScheme(GradientScheme&&) = delete;
    GradientScheme& operator=(GradientScheme&&) = delete;

    /// Destructor
    virtual ~GradientScheme() noexcept = default;

// ****************************** Public Methods ******************************

    /// Calculate gradient at a single cell (scheme-specific)
    [[nodiscard]] virtual Vector cellGradient
    (
        Field field,
        const ScalarField& phi,
        Index cellIdx
    ) const = 0;

    /// Interpolate gradient at a single face
    [[nodiscard]] Vector faceGradient
    (
        Field field,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Index faceIndex
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

// ***************************** Protected Members ****************************

protected:

    /// Construct gradient scheme with mesh context (derived classes only)
    GradientScheme
    (
        const Mesh& mesh,
        const BoundaryConditions& bc
    );

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

    /// Reference to boundary conditions manager
    const BoundaryConditions& bcManager_;

// ****************************** Private Methods *****************************

private:

    /// Distance-weighted linear interpolation of gradients
    Vector averageFaceGradient
    (
        const Face& internalFace,
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

// ****************************** Private Members *****************************

private:

    /// Minimum fraction of ||dPf|| used as normal distance to a boundary face
    /// Prevents gradient amplification beyond ~87 degrees of non-orthogonality
    static constexpr Scalar minNormalFraction_ = S(0.05);
};
