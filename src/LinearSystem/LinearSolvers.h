/******************************************************************************
 * @file LinearSolvers.h
 * @brief Iterative solver for sparse linear systems
 * 
 * @details This class provides configurable solvers with preconditioners for 
 * discretized transport equations. Each LinearSolver instance maintains  
 * independent convergence parameters and preconditioner settings.
 *
 * @class LinearSolver
 *
 * The LinearSolver class provides:
 * - BiCGSTAB solver with Jacobi (diagonal) preconditioner
 * - PCG solver with Incomplete Cholesky for pressure correction
 * - Configurable convergence tolerances
 *****************************************************************************/

#pragma once

#include <string>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/IterativeLinearSolvers>

#include "Scalar.h"


class LinearSolver
{
public:

    /**
     * @brief Construct linear solver with convergence parameters
     * @param fieldName Name of field being solved
     * @param tolerance Absolute residual tolerance
     * @param maxIterations Maximum solver iterations
     */
    LinearSolver
    (
        std::string fieldName,
        Scalar tolerance = S(1e-6),
        int maxIterations = 1000
    );

    /// Copy constructor and assignment - Customized for Eigen members
    LinearSolver(const LinearSolver& other);
    LinearSolver& operator=(const LinearSolver& other);

    /// Move constructor and assignment - deleted: Eigen iterative solver move
    /// is unsound for un-analysed instances (self-referential internal wrapper)
    LinearSolver(LinearSolver&&) = delete;
    LinearSolver& operator=(LinearSolver&&) = delete;

    /// Destructor
    ~LinearSolver() noexcept = default; 

// Preconditioner configuration

    /**
     * @brief Configure Incomplete Cholesky preconditioner for PCG solver
     * @param initialShift Initial shift parameter for numerical stability
     */
    void setICParameters(Scalar initialShift) noexcept
    {
        icInitialShift_ = initialShift;
    }

// Setters

    /**
     * @brief Set absolute residual tolerance
     * @param tol Convergence threshold
     */
    void setTolerance(Scalar tol) noexcept { tolerance_ = tol; }

    /**
     * @brief Set maximum solver iterations
     * @param maxIter Maximum allowed iterations
     */
    void setMaxIterations(int maxIter) noexcept { maxIterations_ = maxIter; }

// Accessors

    /// Get field name for this solver
    [[nodiscard]] const std::string& fieldName() const noexcept
    {
        return fieldName_;
    }

    /// Get absolute tolerance
    [[nodiscard]] Scalar tolerance() const noexcept { return tolerance_; }

    /// Get maximum iterations
    [[nodiscard]] int maxIterations() const noexcept
    {
        return maxIterations_;
    }

    /// Iterations performed by the last solve call (BiCGSTAB or PCG)
    [[nodiscard]] int lastIterations() const noexcept
    {
        return lastIterations_;
    }

    /// Final relative residual reported by the last solve call
    [[nodiscard]] Scalar lastResidual() const noexcept
    {
        return lastResidual_;
    }

// Solver methods

    /**
     * @brief Solve sparse system using BiCGSTAB with Jacobi preconditioner
     * @param x Solution vector
     * @param A Sparse coefficient matrix
     * @param B Right-hand side vector
     *
     * Suitable for non-symmetric systems (momentum, turbulence).
     */
    void solveWithBiCGSTAB
    (
        Eigen::Ref<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x,
        const Eigen::SparseMatrix<Scalar, Eigen::RowMajor>& A,
        const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B
    );

    /**
     * @brief Solve sparse system using PCG with Incomplete Cholesky
     * @param x Solution vector
     * @param A Sparse coefficient matrix (must be SPD)
     * @param B Right-hand side vector
     *
     * Requires symmetric positive definite matrix
     */
    void solveWithPCG
    (
        Eigen::Ref<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x,
        const Eigen::SparseMatrix<Scalar, Eigen::RowMajor>& A,
        const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B
    );
    
private:

    /// Name of field being solved
    std::string fieldName_;

    /// Absolute residual tolerance for convergence
    Scalar tolerance_;

    /// Maximum solver iterations before failure
    int maxIterations_;

    /// Incomplete Cholesky initial shift parameter
    Scalar icInitialShift_ = S(1e-3);

    /// Whether solver parameters and sparsity pattern have been initialized
    bool solverInitialized_ = false;

    /// Iterations performed by the most recent solve call
    int lastIterations_ = 0;

    /// Final residual reported by the most recent solve call
    Scalar lastResidual_ = S(0.0);

    /// BiCGSTAB solver with Jacobi (diagonal) preconditioner — parallel-friendly
    Eigen::BiCGSTAB
    <
        Eigen::SparseMatrix<Scalar, Eigen::RowMajor>,
        Eigen::DiagonalPreconditioner<Scalar>
    >
    bicgstab_;

    /// PCG solver with Incomplete Cholesky preconditioner
    Eigen::ConjugateGradient
    <
        Eigen::SparseMatrix<Scalar, Eigen::RowMajor>,
        Eigen::Lower|Eigen::Upper,
        Eigen::IncompleteCholesky<Scalar>
    >
    pcg_;
};
