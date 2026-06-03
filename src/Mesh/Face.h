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

// ********************************** Headers *********************************

// Standard library headers
#include <optional>
#include <utility>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "OptionalRef.h"
#include "BoundaryPatch.h"
#include "MeshContainers.h"
#include "Integer.h"

// *************************** Forward Declarations ***************************

class Cell;

// *************************** struct FaceIntegrals ***************************

struct FaceIntegrals
{
    Scalar x2 = S(0.0);
    Scalar y2 = S(0.0);
    Scalar z2 = S(0.0);
    Scalar volume = S(0.0);
};

// ******************************** class Face ********************************

class Face
{
public:

    using OptionalIndex = std::optional<Index>;
    using OptionalScalar = std::optional<Scalar>;
    using OptionalVector = std::optional<Vector>;

// ************************* Special Member Functions *************************

    /// Default constructor
    Face() = default;

    /// Constructor for internal faces
    Face
    (
        Index faceIdx,
        IndexList nodes,
        Index owner,
        Index neighbor
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
        Index faceIdx,
        IndexList nodes,
        Index owner
    )
    :
        idx_(faceIdx),
        nodeIndices_(std::move(nodes)),
        ownerCell_(owner),
        neighborCell_(std::nullopt)
    {}

// ****************************** Setter Methods ******************************

    /// Set face identifier
    void setIdx(Index faceIdx) noexcept
    {
        idx_ = faceIdx;
    }

    /// Set owner cell index
    void setOwnerCell(Index owner) noexcept
    {
        ownerCell_ = owner;
    }

    /// Set neighbor cell index
    void setNeighborCell(Index neighbor) noexcept
    {
        neighborCell_ = neighbor;
    }

    /// Set neighbor cell to null
    void setNeighborCell(std::nullopt_t) noexcept
    {
        neighborCell_ = std::nullopt;
    }

    /// Add node index to face connectivity
    void addNodeIndex(Index nodeIdx)
    {
        nodeIndices_.push_back(nodeIdx);
    }

    /// Clear all node indices
    void clearNodeIndices() noexcept
    {
        nodeIndices_.clear();
    }

    /// Set the boundary patch this face belongs to
    void setPatch(const BoundaryPatch& p) noexcept
    {
        patch_ = std::cref(p);
    }

// ***************************** Accessor Methods *****************************

    /// Get face identifier
    [[nodiscard]] Index idx() const noexcept
    {
        return idx_;
    }

    /// Get node connectivity
    [[nodiscard]] IndexListRef nodeIndices() const noexcept
    {
        return nodeIndices_;
    }

    /// Get owner cell index
    [[nodiscard]] Index ownerCell() const noexcept
    {
        return ownerCell_;
    }

    /// Get neighbor cell index
    [[nodiscard]] const OptionalIndex& neighborCell() const noexcept
    {
        return neighborCell_;
    }

    /// Get face centroid
    [[nodiscard]] const Vector& centroid() const noexcept
    {
        return centroid_;
    }

    /// Get face normal vector
    [[nodiscard]] const Vector& normal() const noexcept
    {
        return normal_;
    }

    /// Get face area for flux calculations
    [[nodiscard]] Scalar projectedArea() const noexcept
    {
        return projectedArea_;
    }

    /// Get face contact area
    [[nodiscard]] Scalar contactArea() const noexcept
    {
        return contactArea_;
    }

    /// Get owner cell distance vector
    [[nodiscard]] const Vector& dPf() const noexcept
    {
        return dPf_;
    }

    /// Get neighbor cell distance vector
    [[nodiscard]] const OptionalVector& dNf() const noexcept
    {
        return dNf_;
    }

    /// Get owner cell distance magnitude
    [[nodiscard]] Scalar dPfMag() const noexcept
    {
        return dPfMag_;
    }

    /// Get neighbor cell distance magnitude
    [[nodiscard]] const OptionalScalar& dNfMag() const noexcept
    {
        return dNfMag_;
    }

    /// Get the boundary patch this face belongs to
    [[nodiscard]] const OptionalRef<BoundaryPatch>& patch() const noexcept
    {
        return patch_;
    }

    /// Check if this is a boundary face
    [[nodiscard]] bool isBoundary() const noexcept
    {
        return !neighborCell_.has_value();
    }

// ************************ Geometric Property Methods ************************

    /// Check if geometric properties calculated
    [[nodiscard]] bool geometricPropertiesCalculated() const noexcept
    {
        return geometricPropertiesCalculated_;
    }

    /// Calculate Face centroid, normal, area, and second moment integral
    [[nodiscard]] FaceIntegrals geometricProperties
    (
        NodeListRef allNodes
    );

    /// Check if distance properties calculated
    [[nodiscard]] bool distancesCalculated() const noexcept
    {
        return distancePropertiesCalculated_;
    }

    /// Calculate distance properties of the face
    void distances(CellListRef allCells);

// ****************************** Private Members *****************************

private:

    /// Unique face identifier
    Index idx_ = 0;

    /// Indices of nodes that define this face
    IndexList nodeIndices_;

    /// Index of the owner cell
    Index ownerCell_ = 0;

    /// Index of neighbor cell (nullopt for boundary faces)
    OptionalIndex neighborCell_;

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
    OptionalVector dNf_;

    /// Magnitude of d_Pf
    Scalar dPfMag_ = S(0.0);

    /// Magnitude of d_Nf
    OptionalScalar dNfMag_;

    /// Flag indicating if geometric properties calculated
    bool geometricPropertiesCalculated_ = false;

    /// Flag indicating if distance properties calculated
    bool distancePropertiesCalculated_ = false;

    /// Owning boundary patch (nullopt for internal or unlinked faces)
    OptionalRef<BoundaryPatch> patch_;

// ****************************** Private Methods *****************************

private:

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

// *************************** Non-Member Functions ***************************

/// Stream output operator for Face
std::ostream& operator<<(std::ostream& os, const Face& f);
