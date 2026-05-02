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

#include <vector>
#include <optional>

#include <eigen3/Eigen/SparseCore>

#include "Mesh.h"
#include "BoundaryConditions.h"
#include "TransportEquation.h"


class Matrix
{
public:

    // Eigen type reductions for readability
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;
    using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    /**
     * @brief Constructor for matrix assembly
     * @param mesh Mesh view (nodes, faces, cells)
     * @param boundaryConds Boundary conditions manager
     */
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

    /**
     * @brief Build transport equation matrix
     * @param equation Transport equation data
     */
    void buildMatrix(const TransportEquation& equation);

// Accessor methods

    /**
     * @brief Get assembled sparse matrix A (const)
     * @return Const reference to coefficient matrix
     */
    [[nodiscard]] const SparseMatrix& matrixA() const noexcept
    {
        return matrixA_;
    }

    /**
     * @brief Get assembled sparse matrix A (non-const)
     * @return Mutable reference to coefficient matrix
     */
    [[nodiscard]] SparseMatrix& matrixA() noexcept
    {
        return matrixA_;
    }

    /**
     * @brief Get right-hand side vector b (const)
     * @return Const reference to RHS vector
     */
    [[nodiscard]] const Vec& vectorB() const noexcept
    {
        return vectorB_;
    }

    /**
     * @brief Get right-hand side vector b (non-const)
     * @return Mutable reference to RHS vector
     */
    [[nodiscard]] Vec& vectorB() noexcept
    {
        return vectorB_;
    }

    /**
     * @brief Apply Patankar implicit under-relaxation
     * @param alpha Relaxation factor (0 < alpha <= 1)
     * @param phiPrev Previous iteration field values
     *
     * @details 
     * Modifies the assembled system (Ax=b):
     * - Diagonal: a_P <- a_P / alpha
     * - RHS: b <- b + ((1-alpha)/alpha) * aPOriginal * phiPrev
     */
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

private:

// Private members

    /// Mesh and boundary data references
    const Mesh& mesh_;
    const BoundaryConditions& bcManager_;

    /// Sparse linear system components
    SparseMatrix matrixA_;
    Vec vectorB_;

    /// Triplet storage for efficient sparse matrix assembly
    std::vector<Eigen::Triplet<Scalar>> tripletList_;

    /// Cached face counts for triplet list reservation
    size_t numInternalFaces_ = 0;
    size_t numBoundaryFaces_ = 0;

    /// Relaxation factor from last relax() call (0 = not relaxed)
    Scalar lastRelaxationFactor_ = S(0.0);

    /// Threshold below which f/(1-f) overflows
    Scalar rootSmallValue_ = std::sqrt(smallValue);


// Private methods

    /**
     * @brief Clear matrix and vector for new assembly
     */
    void clear();

    /**
     * @brief Assemble internal face contributions
     * @param face Internal face to process
     * @param equation Transport equation data
     * @param triplets Output triplet list (thread-local buffer)
     * @param localB Output RHS contributions (thread-local buffer)
     */
    void assembleInternalFace
    (
        const Face& face,
        const TransportEquation& equation,
        std::vector<Eigen::Triplet<Scalar>>& triplets,
        Vec& localB
    );

    /**
     * @brief Assemble boundary face contributions
     * @param face Boundary face to process
     * @param equation Transport equation data
     * @param triplets Output triplet list (thread-local buffer)
     * @param localB Output RHS contributions (thread-local buffer)
     */
    void assembleBoundaryFace
    (
        const Face& face,
        const TransportEquation& equation,
        std::vector<Eigen::Triplet<Scalar>>& triplets,
        Vec& localB
    );

    /**
     * @brief Extract scalar boundary value from BC data
     * @param bc Boundary data for the patch/field
     * @param component Component index for vector BCs (0=x, 1=y, 2=z)
     * @return Scalar boundary value
     */
    [[nodiscard]] static Scalar extractBoundaryScalar
    (
        const BoundaryData& bc,
        std::optional<int> component
    ) noexcept;
};

/**
 * @brief Convert a size_t index to Eigen's signed index type
 * @param value Unsigned index from STL/mesh containers
 * @return Equivalent Eigen::Index for sparse-matrix/vector access
 */
[[nodiscard]] inline Eigen::Index eIdx(std::size_t value) noexcept
{
    return static_cast<Eigen::Index>(value);
}
