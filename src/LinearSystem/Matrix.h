/******************************************************************************
 * @file Matrix.h
 * @brief Matrix assembly for finite volume discretization
 *
 * @details This header defines: TransportEquation struct and the Matrix class.
 * TransportEquation bundles all data describing a scalar transport equation
 * (field, convection, diffusion, source, gradients). The Matrix class
 * handles assembly of sparse linear systems (Ax=b) for any transport
 * equation type — momentum, pressure correction, and turbulence.
 *
 * @class Matrix
 *
 * The Matrix class provides:
 * - Unified assembly for all transport equation types
 * - Eigen sparse matrix for efficient storage and solution
 * - Deferred correction for higher-order convection schemes
 * - Non-orthogonal mesh corrections
 * - Implicit under-relaxation (Patankar) for solution stability
 * - Boundary condition integration during assembly
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <vector>
#include <optional>

// External library headers
#include <eigen3/Eigen/SparseCore>

// Project headers
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "TransportEquation.h"

// ******************************* class Matrix *******************************

class Matrix
{
public:

    // Eigen type reductions for readability
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;
    using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

// ************************* Special Member Functions *************************

    /// Constructor for matrix assembly
    Matrix
    (
        const Mesh& mesh,
        const BoundaryConditions& boundaryConds
    ) noexcept;

    /// Copy constructor and assignment - Not copyable (const T& members)
    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    Matrix(Matrix&&) = delete;
    Matrix& operator=(Matrix&&) = delete;

    /// Destructor
    ~Matrix() noexcept = default;

// ****************************** Public Methods ******************************

    /// Build transport equation matrix
    void buildMatrix(const TransportEquation& equation);

// ***************************** Accessor Methods *****************************

    /// Get assembled sparse matrix A (const)
    [[nodiscard]] const SparseMatrix& matrixA() const noexcept
    {
        return matrixA_;
    }

    /// Get assembled sparse matrix A (non-const)
    [[nodiscard]] SparseMatrix& matrixA() noexcept
    {
        return matrixA_;
    }

    /// Get right-hand side vector b (const)
    [[nodiscard]] const Vec& vectorB() const noexcept
    {
        return vectorB_;
    }

    /// Get right-hand side vector b
    [[nodiscard]] Vec& vectorB() noexcept
    {
        return vectorB_;
    }

    /// Apply Patankar implicit under-relaxation
    void relax(Scalar alpha, const ScalarField& phiPrev);

    /**
     * @brief Fix matrix rows to impose known cell values
     *
     * @details
     * Replaces the equation for each constrained cell with
     * diag * phi[i] = diag * value[i], and transfers the
     * known-value coupling to unconstrained neighbors' RHS.
     * Call after relax() but before solve().
     *
     * @param cellIndices Indices of cells to constrain
     * @param values Prescribed values for those cells
     */
    void setValues
    (
        std::span<const size_t> cellIndices,
        std::span<const Scalar> values,
        std::span<const Scalar> fractions = {}
    );

// ****************************** Private Members *****************************

private:

    /// Mesh data reference
    const Mesh& mesh_;

    /// Boundary data references
    const BoundaryConditions& bcManager_;

    /// Sparse linear system components
    SparseMatrix matrixA_;
    Vec vectorB_;

    /// Triplet storage for efficient sparse matrix assembly
    std::vector<Eigen::Triplet<Scalar>> tripletList_;

    /// Per-thread triplet lists for parallel face-assembly
    std::vector<std::vector<Eigen::Triplet<Scalar>>> perThreadTriplets_;

    /// Per-thread RHS contributions for parallel face-assembly scatter.
    std::vector<Vec> perThreadB_;

    /// Cached face counts for triplet list reservation
    size_t numInternalFaces_ = 0;
    size_t numBoundaryFaces_ = 0;

    /// Relaxation factor from last relax() call (0 = not relaxed)
    Scalar lastRelaxationFactor_ = S(0.0);

    /// Threshold below which f/(1-f) overflows
    inline static const Scalar rootSmallValue_ = std::sqrt(smallValue);


// ****************************** Private Methods *****************************

private:

    /// Clear matrix and vector for new assembly
    void clear();

    /// Assemble internal face contributions
    void assembleInternalFace
    (
        const Face& face,
        const TransportEquation& equation,
        std::vector<Eigen::Triplet<Scalar>>& triplets,
        Vec& localB
    ) const;

    /// Assemble boundary face contributions
    void assembleBoundaryFace
    (
        const Face& face,
        const TransportEquation& equation,
        std::vector<Eigen::Triplet<Scalar>>& triplets,
        Vec& localB
    ) const;

    /// Resolve the diffusion coefficient at a face
    [[nodiscard]] Scalar faceDiffusionCoefficient
    (
        const Face& face,
        const TransportEquation& equation
    ) const;
};

// *************************** Non-Member Functions ***************************

/// Convert a size_t index to Eigen's signed index type
[[nodiscard]] inline Eigen::Index eIdx(std::size_t value) noexcept
{
    return static_cast<Eigen::Index>(value);
}
