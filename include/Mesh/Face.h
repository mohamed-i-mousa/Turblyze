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
    Scalar x2 = 0.0;
    Scalar y2 = 0.0;
    Scalar z2 = 0.0;
    Scalar volume = 0.0;
};


class Face
{
public:

    /// Default constructor
    Face() = default;

    /**
     * @brief Constructor for internal faces
     * @param faceIdx Unique face identifier
     * @param nodes Indices of face nodes
     * @param owner Index of owner cell
     * @param neighbor Index of neighbor cell
     */
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

    /**
     * @brief Constructor for boundary faces
     * @param faceIdx Unique face identifier
     * @param nodes Indices of face nodes
     * @param owner Index of owner cell
     */
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

    /**
     * @brief Set face identifier
     * @param faceIdx Unique face ID
     */
    void setIdx(size_t faceIdx) noexcept { idx_ = faceIdx; }

    /**
     * @brief Set owner cell index
     * @param owner Index of owner cell
     */
    void setOwnerCell(size_t owner) noexcept { ownerCell_ = owner; }

    /**
     * @brief Set neighbor cell index
     * @param neighbor Index of neighbor cell
     */
    void setNeighborCell(size_t neighbor) noexcept
    {
        neighborCell_ = neighbor;
    }

    /**
     * @brief Set neighbor cell to null (boundary face)
     */
    void setNeighborCell(std::nullopt_t) noexcept
    {
        neighborCell_ = std::nullopt;
    }

    /**
     * @brief Add node index to face connectivity
     * @param nodeIdx Index of node to add
     */
    void addNodeIndex(size_t nodeIdx) { nodeIndices_.push_back(nodeIdx); }

    /**
     * @brief Clear all node indices
     */
    void clearNodeIndices() noexcept { nodeIndices_.clear(); }

    /**
     * @brief Set the boundary patch this face belongs to
     * @param p Reference to the owning boundary patch
     */
    void setPatch(const BoundaryPatch& p) noexcept { patch_ = std::cref(p); }

// Accessor methods

    /**
     * @brief Get face identifier
     * @return Unique face ID
     */
    [[nodiscard]] size_t idx() const noexcept { return idx_; }

    /**
     * @brief Get node connectivity
     * @return Vector of node indices
     */
    [[nodiscard]] std::span<const size_t> nodeIndices() const noexcept
    {
        return nodeIndices_;
    }

    /**
     * @brief Get owner cell index
     * @return Index of owner cell
     */
    [[nodiscard]] size_t ownerCell() const noexcept { return ownerCell_; }

    /**
     * @brief Get neighbor cell index
     * @return Optional neighbor cell index
     */
    [[nodiscard]] const std::optional<size_t>& neighborCell() const noexcept
    {
        return neighborCell_;
    }

    /**
     * @brief Get face centroid
     * @return Face center coordinates
     */
    [[nodiscard]] const Vector& centroid() const noexcept { return centroid_; }

    /**
     * @brief Get face normal vector
     * @return Unit normal vector
     */
    [[nodiscard]] const Vector& normal() const noexcept { return normal_; }

    /**
     * @brief Get face area (projected area) for flux calculations
     * @return Face area magnitude (projected area for non-planar faces)
     */
    [[nodiscard]] Scalar projectedArea() const noexcept
    {
        return projectedArea_;
    }

    /**
     * @brief Get face contact area (actual wetted surface area)
     * @return Contact area (sum of sub-triangle areas for non-planar faces)
     * @note Used for wall shear stress, heat transfer, and friction drag
     */
    [[nodiscard]] Scalar contactArea() const noexcept { return contactArea_; }

    /**
     * @brief Get owner cell distance vector
     * @return Vector from owner to face
     */
    [[nodiscard]] const Vector& dPf() const noexcept { return dPf_; }

    /**
     * @brief Get neighbor cell distance vector
     * @return Optional vector from neighbor to face
     */
    [[nodiscard]] const std::optional<Vector>& dNf() const noexcept
    {
        return dNf_;
    }

    /**
     * @brief Get owner cell distance magnitude
     * @return Distance from owner to face
     */
    [[nodiscard]] Scalar dPfMag() const noexcept { return dPfMag_; }

    /**
     * @brief Get neighbor cell distance magnitude
     * @return Optional distance from neighbor to face
     */
    [[nodiscard]] const std::optional<Scalar>& dNfMag() const noexcept
    {
        return dNfMag_;
    }

    /**
     * @brief Check if geometric properties calculated
     * @return True if geometry computed
     */
    [[nodiscard]] bool geometricPropertiesCalculated() const noexcept
    {
        return geometricPropertiesCalculated_;
    }

    /**
     * @brief Check if distance properties calculated
     * @return True if distance vectors and magnitudes computed
     */
    [[nodiscard]] bool distancePropertiesCalculated() const noexcept
    {
        return distancePropertiesCalculated_;
    }

    /**
     * @brief Calculate geometric properties of the face
     *
     * @details
     * - Calculates face area, centroid, normal vector, and second moment
     *   integrals.
     * - For triangles, uses direct cross product. For polygons, decomposes
     *   into triangles.
     * - Sets geometricPropertiesCalculated flag when success.
     *
     * @param allNodes Vector of all mesh nodes
     * @note Terminates the program if node index is invalid
     * @note Terminates the program if face has zero area (collinear nodes)
     * @return FaceIntegrals for cell volume/centroid computation
     */
    [[nodiscard]] FaceIntegrals calculateGeometricProperties
    (
        std::span<const Vector> allNodes
    );

    /**
     * @brief Calculate distance properties of the face
     *
     * @details
     * - Calculates distance vectors, magnitudes, and unit vectors
     *   from cell centers to face center. For boundary faces,
     *   only owner cell distances are calculated.
     *
     * @param cellCentroids Centroid of every cell indexed by cell index
     */
    void calculateDistanceProperties(std::span<const Vector> cellCentroids);

    /**
     * @brief Check if this is a boundary face
     * @return True if face is on domain boundary
     */
    [[nodiscard]] bool isBoundary() const noexcept
    {
        return !neighborCell_.has_value();
    }

    /**
     * @brief Get the boundary patch this face belongs to
     * @return Optional reference to owning patch (nullopt if unlinked)
     */
    [[nodiscard]] const OptionalRef<BoundaryPatch>& patch() const noexcept
    {
        return patch_;
    }

    /**
     * @brief Flip the face normal direction
     */
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
    Scalar projectedArea_ = 0.0;

    /// Contact area (For shear stress calculations)
    Scalar contactArea_ = 0.0;

    /// Distance vector from owner cell center to face center
    Vector dPf_;

    /// Distance vector from neighbor cell center to face center
    std::optional<Vector> dNf_;

    /// Magnitude of d_Pf
    Scalar dPfMag_ = 0.0;

    /// Magnitude of d_Nf
    std::optional<Scalar> dNfMag_;

    /// Flag indicating if geometric properties calculated
    bool geometricPropertiesCalculated_ = false;

    /// Flag indicating if distance properties calculated
    bool distancePropertiesCalculated_ = false;

    /// Owning boundary patch (nullopt for internal or unlinked faces)
    OptionalRef<BoundaryPatch> patch_;

    /**
     * @brief Symmetric second-moment polynomial for triangle integration
     *
     * @details Evaluates a² + b² + c² + ab + ac + bc, which appears in the
     * divergence-theorem integral of x² over a triangle whose vertices have
     * coordinate values a, b, c along a given axis:
     * ∫∫_triangle x² dA = (area / 6) × secondMoment(x₁, x₂, x₃)
     *
     * @param a First vertex coordinate along one axis
     * @param b Second vertex coordinate along one axis
     * @param c Third vertex coordinate along one axis
     * @return a² + b² + c² + ab + ac + bc
     */
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

/**
 * @brief Stream output operator for Face
 * @param os Output stream
 * @param f Face to output
 * @return Reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Face& f);
