#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <vector>
#include <iostream>
#include <cmath>
#include "Scalar.h"

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
 * @brief Solve symmetric linear system Ax = b using Preconditioned Conjugate Gradient method
 * @param x Solution vector (input: initial guess, output: solution)
 * @param A Symmetric coefficient matrix
 * @param B Right-hand side vector
 * @param tolerance Convergence tolerance
 * @param maxIterations Maximum number of iterations
 * @param fieldName Name of the field being solved (for logging)
 * @return true if converged, false otherwise
 */
bool PCG
(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B,
    Scalar tolerance,
    int maxIterations,
    const std::string& fieldName
);

}

#endif