/******************************************************************************
 * @file Matrix.hpp
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

#include "Cell.hpp"
#include "Face.hpp"
#include "BoundaryConditions.hpp"
#include "TransportEquation.hpp"


class Matrix
{
public:

    /**
     * @brief Constructor for matrix assembly
     * @param faces Reference to face data
     * @param cells Reference to cell data
     * @param boundaryConds Boundary conditions manager
     */
    Matrix
    (
        const std::span<const Face> faces,
        const std::span<const Cell> cells,
        const BoundaryConditions& boundaryConds
    ) noexcept;

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
    const Eigen::SparseMatrix<Scalar>& matrixA() const noexcept
    {
        return matrixA_;
    }

    /**
     * @brief Get assembled sparse matrix A (non-const)
     * @return Mutable reference to coefficient matrix
     */
    Eigen::SparseMatrix<Scalar>& matrixA() noexcept
    {
        return matrixA_;
    }

    /**
     * @brief Get right-hand side vector b (const)
     * @return Const reference to RHS vector
     */
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&
    vectorB() const noexcept
    {
        return vectorB_;
    }

    /**
     * @brief Get right-hand side vector b (non-const)
     * @return Mutable reference to RHS vector
     */
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& vectorB() noexcept
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
    std::span<const Face> allFaces_;
    std::span<const Cell> allCells_;
    const BoundaryConditions& bcManager_;

    /// Sparse linear system components
    Eigen::SparseMatrix<Scalar> matrixA_;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vectorB_;

    /// Triplet storage for efficient sparse matrix assembly
    std::vector<Eigen::Triplet<Scalar>> tripletList_;

    /// Cached face counts for triplet list reservation
    size_t numInternalFaces_;
    size_t numBoundaryFaces_;

    /// Relaxation factor from last relax() call (0 = not relaxed)
    Scalar lastRelaxationFactor_ = S(0.0);

// Private methods

    /**
     * @brief Clear matrix and vector for new assembly
     */
    void clear();

    /**
     * @brief Reserve triplet list capacity based on mesh topology
     */
    void reserveTripletList();

    /**
     * @brief Assemble internal face contributions
     * @param face Internal face to process
     * @param equation Transport equation data
     */
    void assembleInternalFace
    (
        const Face& face,
        const TransportEquation& equation
    );

    /**
     * @brief Assemble boundary face contributions
     * @param face Boundary face to process
     * @param equation Transport equation data
     */
    void assembleBoundaryFace
    (
        const Face& face,
        const TransportEquation& equation
    );

    /**
     * @brief Extract scalar boundary value from BC data
     * @param bc Boundary data for the patch/field
     * @param component Component index for vector BCs (0=x, 1=y, 2=z)
     * @return Scalar boundary value
     */
    static Scalar extractBoundaryScalar
    (
        const BoundaryData& bc,
        std::optional<int> component
    ) noexcept;
};
