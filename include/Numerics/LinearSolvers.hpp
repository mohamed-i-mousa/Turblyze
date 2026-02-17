/******************************************************************************
 * @file LinearSolvers.hpp
 * @brief Iterative solver for sparse linear systems
 * 
 * This class provides configurable solvers with preconditioners for 
 * discretized transport equations. Each LinearSolver instance maintains  
 * independent convergence parameters and preconditioner settings.
 *
 * @class LinearSolver
 *
 * The LinearSolver class provides:
 * - BiCGSTAB solver with ILUT preconditioner for momentum and turbulence
 * - PCG solver with Incomplete Cholesky for pressure correction
 * - Configurable convergence tolerances and preconditioner parameters
 * - Optional exact residual reporting with convergence diagnostics
 *****************************************************************************/

#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <string>
#include <iostream>
#include <cmath>
#include <optional>

#include "Scalar.hpp"

class LinearSolver
{
public:

    /**
     * @brief Construct linear solver with convergence parameters
     * @param fieldName Name of field being solved
     * @param tolerance Absolute residual tolerance
     * @param maxIterations Maximum solver iterations
     * @param relTol Relative residual tolerance (final/initial)
     * @param verbose Enable detailed iteration output
     */
    LinearSolver
    (
        const std::string& fieldName,
        Scalar tolerance = S(1e-6),
        int maxIterations = 1000,
        Scalar relTol = S(0.0),
        bool verbose = false
    );

// Preconditioner configuration

    /**
     * @brief Configure ILUT preconditioner for BiCGSTAB solver
     * @param fillFactor Maximum fill-in factor
     * @param dropTol Drop tolerance for small entries during factorization
     */
    void setILUTParameters(int fillFactor, Scalar dropTol)
    {
        ilutFillFactor_ = fillFactor;
        ilutDropTol_ = dropTol;
    }

    /**
     * @brief Configure Incomplete Cholesky preconditioner for PCG solver
     * @param initialShift Initial shift parameter for numerical stability
     */
    void setICParameters(Scalar initialShift)
    {
        icInitialShift_ = initialShift;
    }

// Setters

    /**
     * @brief Set absolute residual tolerance
     * @param tol Convergence threshold
     */
    void setTolerance(Scalar tol) { tolerance_ = tol; }

    /**
     * @brief Set maximum solver iterations
     * @param maxIter Maximum allowed iterations
     */
    void setMaxIterations(int maxIter) { maxIterations_ = maxIter; }

    /**
     * @brief Set relative residual tolerance
     * @param relTol Ratio of final to initial residual for convergence
     */
    void setRelTol(Scalar relTol) { relTol_ = relTol; }

    /**
     * @brief Enable or disable verbose output
     * @param verbose True to print iteration statistics
     */
    void setVerbose(bool verbose) { verbose_ = verbose; }

    /**
     * @brief Enable or disable exact residual computation
     * @param compute True to compute B-A*x residuals each solve
     */
    void setComputeResiduals(bool compute)
    {
        computeResiduals_ = compute;
    }

// Accessors

    /// Get field name for this solver
    const std::string& fieldName() const { return fieldName_; }

    /// Get absolute tolerance
    Scalar tolerance() const { return tolerance_; }

    /// Get maximum iterations
    int maxIterations() const { return maxIterations_; }

    /// Get relative tolerance
    Scalar relTol() const { return relTol_; }

    /// Check if exact residual computation is enabled
    bool computeResiduals() const
    {
        return computeResiduals_;
    }

// Solver methods

    /**
     * @brief Solve sparse system using BiCGSTAB with ILUT
     * @param x Solution vector
     * @param A Sparse coefficient matrix
     * @param B Right-hand side vector
     * @return True unless catastrophic failure detected
     *
     * Suitable for non-symmetric systems (momentum, turbulence).
     */
    bool solveWithBiCGSTAB
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
     * @return True unless catastrophic failure detected
     *
     * Requires symmetric positive definite matrix
     */
    bool solveWithPCG
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

    /// Relative residual tolerance (final / initial)
    Scalar relTol_;

    /// Enable detailed iteration output
    bool verbose_;

    /// ILUT fill factor (controls memory vs accuracy)
    int ilutFillFactor_ = 5;

    /// ILUT drop tolerance for small entries
    Scalar ilutDropTol_ = S(1e-6);

    /// Incomplete Cholesky initial shift parameter
    Scalar icInitialShift_ = S(1e-2);

    /// Enable exact B-A*x residual computation
    bool computeResiduals_ = true;

// Shared helpers

    /**
     * @brief Container for residual metrics during solve
     */
    struct ResidualInfo
    {
        Scalar initialResidual;        ///< Initial RMS residual
        Scalar finalResidual;          ///< Final RMS residual
        Scalar residualRatio;          ///< finalResidual / initialResidual
    };

    /**
     * @brief Compute RMS residual for convergence check
     * @param residualVector Residual vector (Ax - b)
     * @return Root mean square of residuals
     */
    Scalar rmsResidual
    (
        const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&
            residualVector
    ) const;

    /**
     * @brief Convert our tolerance to Eigen's convention
     * @param initialResidual RMS of initial residual
     * @param bNorm L2 norm of RHS vector
     * @param n Number of unknowns
     * @return Tolerance for Eigen's setTolerance()
     *
     * Maps max(absTol, relTol * initRMS) to Eigen's
     * |r|/|b| criterion via sqrt(n) * effectiveTol / |b|.
     */
    Scalar computeEigenTolerance
    (
        Scalar initialResidual,
        Scalar bNorm,
        Eigen::Index n
    ) const;

    /**
     * @brief Report solver convergence diagnostics
     * @param iterations Number of iterations performed
     * @param eigenInfo Eigen's internal convergence status
     * @param residuals Optional residual metrics
     *
     * Reports iteration count and Eigen status. When residuals
     * are provided, also reports initial/final RMS residuals,
     * residual ratio vs relTol, and final residual vs tolerance.
     */
    void reportConvergence
    (
        int iterations,
        Eigen::ComputationInfo eigenInfo,
        const std::optional<ResidualInfo>& residuals
    ) const;
};

#endif // LINEAR_SOLVER_HPP