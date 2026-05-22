/******************************************************************************
 * @file Face.h
 * @brief Represents a face in the computational mesh
 *
 * @details This header defines the Face class, which is fundamental in
 * the finite volume discretization.
 * A face represents a surface defined by a sequence of nodes (vertices)
 * and serves as the boundary between two control volumes (cells) or between a
 * cell and the domain boundary.
 *
 * @struct FaceIntegrals
 * - Stores second moment integrals and volume contribution from a face
 *   used in cell geometric property calculations.
 *
 * @class Face
 * - Connectivity (nodes, owner cell, neighbor cell)
 * - Face properties (centroid, area, normal vector)
 * - Distance vectors for interpolations and gradient calculations
 * - Boundary handling (internal and boundary faces)
 *****************************************************************************/

#pragma once

#include <vector>
#include <string>
#include <optional>
#include <span>

#include "Scalar.h"
#include "Vector.h"
#include "OptionalRef.h"
#include "BoundaryPatch.h"


struct FaceIntegrals
{
    Scalar x2 = S(0.0);
    Scalar y2 = S(0.0);
    Scalar z2 = S(0.0);
    Scalar volume = S(0.0);
};


class Face
{
public:

    /// Default constructor
    Face() = default;

    /// Constructor for internal faces
    Face
    (
        size_t faceIdx,
        std::vector<size_t> nodes,
        size_t owner,
        size_t neighbor
    )
    :
        idx_(faceIdx),
        nodeIndices_(std::move(nodes)),
        ownerCell_(owner),
        neighborCell_(neighbor)
    {}

    /// Constructor for boundary faces
    Face
    (
        size_t faceIdx,
        std::vector<size_t> nodes,
        size_t owner
    )
    :
        idx_(faceIdx),
        nodeIndices_(std::move(nodes)),
        ownerCell_(owner),
        neighborCell_(std::nullopt)
    {}

// Setter methods

    /// Set face identifier
    void setIdx(size_t faceIdx) noexcept { idx_ = faceIdx; }

    /// Set owner cell index
    void setOwnerCell(size_t owner) noexcept { ownerCell_ = owner; }

    /// Set neighbor cell index
    void setNeighborCell(size_t neighbor) noexcept
    {
        neighborCell_ = neighbor;
    }

    /// Set neighbor cell to null
    void setNeighborCell(std::nullopt_t) noexcept
    {
        neighborCell_ = std::nullopt;
    }

    /// Add node index to face connectivity
    void addNodeIndex(size_t nodeIdx) { nodeIndices_.push_back(nodeIdx); }

    /// Clear all node indices
    void clearNodeIndices() noexcept { nodeIndices_.clear(); }

    /// Set the boundary patch this face belongs to
    void setPatch(const BoundaryPatch& p) noexcept { patch_ = std::cref(p); }

// Accessor methods

    /// Get face identifier
    [[nodiscard]] size_t idx() const noexcept { return idx_; }

    /// Get node connectivity
    [[nodiscard]] std::span<const size_t> nodeIndices() const noexcept
    {
        return nodeIndices_;
    }

    /// Get owner cell index
    [[nodiscard]] size_t ownerCell() const noexcept { return ownerCell_; }

    /// Get neighbor cell index
    [[nodiscard]] const std::optional<size_t>& neighborCell() const noexcept
    {
        return neighborCell_;
    }

    /// Get face centroid
    [[nodiscard]] const Vector& centroid() const noexcept { return centroid_; }

    /// Get face normal vector
    [[nodiscard]] const Vector& normal() const noexcept { return normal_; }

    /// Get face area for flux calculations
    [[nodiscard]] Scalar projectedArea() const noexcept
    {
        return projectedArea_;
    }

    /// Get face contact area
    [[nodiscard]] Scalar contactArea() const noexcept { return contactArea_; }

    /// Get owner cell distance vector
    [[nodiscard]] const Vector& dPf() const noexcept { return dPf_; }

    /// Get neighbor cell distance vector
    [[nodiscard]] const std::optional<Vector>& dNf() const noexcept
    {
        return dNf_;
    }

    /// Get owner cell distance magnitude
    [[nodiscard]] Scalar dPfMag() const noexcept { return dPfMag_; }

    /// Get neighbor cell distance magnitude
    [[nodiscard]] const std::optional<Scalar>& dNfMag() const noexcept
    {
        return dNfMag_;
    }

    /// Check if geometric properties calculated
    [[nodiscard]] bool geometricPropertiesCalculated() const noexcept
    {
        return geometricPropertiesCalculated_;
    }

    /// Check if distance properties calculated
    [[nodiscard]] bool distancePropertiesCalculated() const noexcept
    {
        return distancePropertiesCalculated_;
    }

    /// Calculate Face centroid, normal, area, and second moment integral
    [[nodiscard]] FaceIntegrals calculateGeometricProperties
    (
        std::span<const Vector> allNodes
    );

    /// Calculate distance properties of the face
    void calculateDistanceProperties(std::span<const Vector> cellCentroids);

    /// Check if this is a boundary face
    [[nodiscard]] bool isBoundary() const noexcept
    {
        return !neighborCell_.has_value();
    }

    /// Get the boundary patch this face belongs to
    [[nodiscard]] const OptionalRef<BoundaryPatch>& patch() const noexcept
    {
        return patch_;
    }

    /// Flip the face normal direction
    void flipNormal() noexcept { normal_ *= S(-1.0); }

private:

    /// Unique face identifier
    size_t idx_ = 0;

    /// Indices of nodes that define this face
    std::vector<size_t> nodeIndices_;

    /// Index of the owner cell
    size_t ownerCell_ = 0;

    /// Index of neighbor cell (nullopt for boundary faces)
    std::optional<size_t> neighborCell_;

    /// Face geometric centroid
    Vector centroid_;

    /// Face normal vector (unit vector)
    Vector normal_;

    /// Face area (projected area for flux calculations)
    Scalar projectedArea_ = S(0.0);

    /// Contact area (For shear stress calculations)
    Scalar contactArea_ = S(0.0);

    /// Distance vector from owner cell center to face center
    Vector dPf_;

    /// Distance vector from neighbor cell center to face center
    std::optional<Vector> dNf_;

    /// Magnitude of d_Pf
    Scalar dPfMag_ = S(0.0);

    /// Magnitude of d_Nf
    std::optional<Scalar> dNfMag_;

    /// Flag indicating if geometric properties calculated
    bool geometricPropertiesCalculated_ = false;

    /// Flag indicating if distance properties calculated
    bool distancePropertiesCalculated_ = false;

    /// Owning boundary patch (nullopt for internal or unlinked faces)
    OptionalRef<BoundaryPatch> patch_;

    /// Symmetric second-moment polynomial for triangle integration
    /// Evaluates a² + b² + c² + ab + ac + bc
    /// ∫∫_triangle x² dA = (area / 6) × secondMoment(x₁, x₂, x₃)
    [[nodiscard]] static Scalar secondMoment
    (
        Scalar a,
        Scalar b,
        Scalar c
    ) noexcept
    {
        return a*a + b*b + c*c + a*b + a*c + b*c;
    }

};

/// Stream output operator for Face
std::ostream& operator<<(std::ostream& os, const Face& f);
