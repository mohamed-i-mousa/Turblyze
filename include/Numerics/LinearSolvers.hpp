/******************************************************************************
 * @file LinearSolvers.h
 * @brief Linear system solvers for sparse matrices in CFD applications
 * 
 * This header provides robust iterative linear solvers optimized for sparse
 * matrices arising from finite volume discretization of transport equations.
 * The solvers are built on the Eigen library and include preconditioned
 * iterative methods with convergence monitoring and adaptive tolerance
 * adjustment for numerical stability.
 * 
 * @namespace LinearSolvers
 * 
 * The LinearSolvers namespace provides:
 * - BiCGSTAB solver with incomplete LU preconditioning
 * - Conjugate gradient solver for symmetric positive definite systems
 * - Adaptive tolerance and maximum iteration controls
 * - Convergence diagnostics and residual monitoring
 * 
 * Key features:
 * - Eigen-based implementation for optimal performance
 * - ILUT preconditioning for faster convergence
 * - Configurable precision (single/double) via Scalar type
 * - Robust error handling and convergence detection
 *****************************************************************************/

#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <vector>
#include <iostream>
#include <cmath>
#include "Scalar.hpp"

/**
 * @brief Namespace containing linear system solvers
 */
namespace LinearSolvers 
{

/**
 * @brief Solve linear system Ax = b using BiCGSTAB method
 * @param x Solution vector (input: initial guess, output: solution)
 * @param A Coefficient matrix
 * @param B Right-hand side vector
 * @param tolerance Convergence tolerance
 * @param maxIterations Maximum number of iterations
 * @param fieldName Name of the field being solved (for logging)
 * @return true if converged, false otherwise
 */
bool BiCGSTAB
(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B,
    Scalar tolerance,
    int maxIterations,
    const std::string& fieldName
);

/**
 * @brief Solve symmetric linear system using Preconditioned CG method
 * @param x Solution vector (input: initial guess, output: solution)
 * @param A Symmetric coefficient matrix
 * @param B Right-hand side vector
 * @param tolerance Absolute convergence tolerance
 * @param maxIterations Maximum number of iterations
 * @param fieldName Name of the field being solved (for logging)
 * @param relTol Relative tolerance (convergence if residual reduces by this factor)
 * @return true if converged, false otherwise
 */
bool PCG
(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B,
    Scalar tolerance,
    int maxIterations,
    const std::string& fieldName,
    Scalar relTol = S(0.05)
);

}

#endif