/******************************************************************************
 * @file LinearSolvers.hpp
 * @brief Iterative solver for sparse linear systems
 * 
 * @details This class provides configurable solvers with preconditioners for 
 * discretized transport equations. Each LinearSolver instance maintains  
 * independent convergence parameters and preconditioner settings.
 *
 * @class LinearSolver
 *
 * The LinearSolver class provides:
 * - BiCGSTAB solver with ILUT preconditioner for momentum and turbulence
 * - PCG solver with Incomplete Cholesky for pressure correction
 * - Configurable convergence tolerances and preconditioner parameters
 *****************************************************************************/

#pragma once

#include <string>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/IterativeLinearSolvers>

#include "Scalar.hpp"

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

    /**
     * @brief Copy configuration to a new solver instance
     * @param other Source solver to copy configuration from
     * @note Copies tolerance, iterations, preconditioner parameters, and
     * debug flag only. Eigen solver state (bicgstab_, pcg_) is not copied —
     * the new instance builds its own on first solve call.
     */
    LinearSolver(const LinearSolver& other);
    LinearSolver& operator=(const LinearSolver& other);

    /// Deleted — Eigen types are not moveable
    LinearSolver(LinearSolver&&) = delete;
    LinearSolver& operator=(LinearSolver&&) = delete;

// Preconditioner configuration

    /**
     * @brief Configure ILUT preconditioner for BiCGSTAB solver
     * @param fillFactor Maximum fill-in factor
     * @param dropTol Drop tolerance for small entries during factorization
     */
    void setILUTParameters(int fillFactor, Scalar dropTol) noexcept
    {
        ilutFillFactor_ = fillFactor;
        ilutDropTol_ = dropTol;
    }

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

    /**
     * @brief Enable or disable debug logging after each solve
     * @param debug True to print iteration count and residual
     */
    void setDebug(bool debug) noexcept { debug_ = debug; }

// Accessors

    /// Get field name for this solver
    const std::string& fieldName() const noexcept
    {
        return fieldName_;
    }

    /// Get absolute tolerance
    Scalar tolerance() const noexcept { return tolerance_; }

    /// Get maximum iterations
    int maxIterations() const noexcept
    {
        return maxIterations_;
    }

// Solver methods

    /**
     * @brief Solve sparse system using BiCGSTAB with ILUT
     * @param x Solution vector
     * @param A Sparse coefficient matrix
     * @param B Right-hand side vector
     *
     * Suitable for non-symmetric systems (momentum, turbulence).
     */
    void solveWithBiCGSTAB
    (
        Eigen::Ref<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x,
        const Eigen::SparseMatrix<Scalar>& A,
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
        const Eigen::SparseMatrix<Scalar>& A,
        const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B
    );
    
private:

    /// Name of field being solved
    std::string fieldName_;

    /// Absolute residual tolerance for convergence
    Scalar tolerance_;

    /// Maximum solver iterations before failure
    int maxIterations_;

    /// ILUT fill factor (controls memory vs accuracy)
    int ilutFillFactor_ = 5;

    /// ILUT drop tolerance for small entries
    Scalar ilutDropTol_ = S(1e-3);

    /// Incomplete Cholesky initial shift parameter
    Scalar icInitialShift_ = S(1e-3);

    /// Whether solver parameters and sparsity pattern have been initialized
    bool solverInitialized_ = false;

    /// Enable debug logging of iteration count and residual after each solve
    bool debug_ = false;

    /// BiCGSTAB solver with ILUT preconditioner
    Eigen::BiCGSTAB
    <Eigen::SparseMatrix<Scalar>, Eigen::IncompleteLUT<Scalar>>
    bicgstab_;

    /// PCG solver with Incomplete Cholesky preconditioner
    Eigen::ConjugateGradient
    <
        Eigen::SparseMatrix<Scalar>,
        Eigen::Lower|Eigen::Upper,
        Eigen::IncompleteCholesky<Scalar>
    >
    pcg_;
};
