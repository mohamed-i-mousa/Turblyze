/******************************************************************************
 * @file LinearSolvers.cpp
 * @brief LinearSolver class implementation
 *****************************************************************************/

#include "LinearSolvers.hpp"

#include <iostream>

#include "ErrorHandler.hpp"


// ************************* Constructor & Setters ************************

LinearSolver::LinearSolver
(
    std::string fieldName,
    Scalar tolerance,
    int maxIterations
)
:
    fieldName_(std::move(fieldName)),
    tolerance_(tolerance),
    maxIterations_(maxIterations)
{}


// *********************************** Copy ***********************************

LinearSolver::LinearSolver(const LinearSolver& other)
:   
    fieldName_(other.fieldName_),
    tolerance_(other.tolerance_),
    maxIterations_(other.maxIterations_),
    ilutFillFactor_(other.ilutFillFactor_),
    ilutDropTol_(other.ilutDropTol_),
    icInitialShift_(other.icInitialShift_),
    debug_(other.debug_)
{}

LinearSolver& LinearSolver::operator=(const LinearSolver& other)
{
    if (this != &other)
    {
        fieldName_       = other.fieldName_;
        tolerance_       = other.tolerance_;
        maxIterations_   = other.maxIterations_;
        ilutFillFactor_  = other.ilutFillFactor_;
        ilutDropTol_     = other.ilutDropTol_;
        icInitialShift_  = other.icInitialShift_;
        debug_           = other.debug_;
        solverInitialized_ = false;
    }
    return *this;
}


// ***************************** BiCGSTAB Solver *****************************

void LinearSolver::solveWithBiCGSTAB
(
    Eigen::Ref<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B
)
{
    if (x.size() != B.size())
    {
        FatalError("BiCGSTAB: x and B size mismatch");
    }

    bicgstab_.setMaxIterations(maxIterations_);
    bicgstab_.setTolerance(tolerance_);
    bicgstab_.preconditioner().setFillfactor(ilutFillFactor_);
    bicgstab_.preconditioner().setDroptol(ilutDropTol_);

    // Sparsity analysis — run only once
    if (!solverInitialized_)
    {
        bicgstab_.analyzePattern(A);
        solverInitialized_ = true;
    }

    // Numeric factorization — coefficients change
    bicgstab_.factorize(A);

    if (bicgstab_.info() != Eigen::Success)
    {
        Warning
        (
            "Field '" + fieldName_
          + "': BiCGSTAB factorization failed!"
        );
        return;
    }

    x = bicgstab_.solveWithGuess(B, x);

    if (debug_)
    {
        std::cout
            << "  [" << fieldName_ << "] BiCGSTAB: "
            << bicgstab_.iterations() << " iterations"
            << ", residual = " << std::scientific
            << bicgstab_.error() << std::fixed
            << std::endl;
    }
}

// ******************************* PCG Solver *******************************

void LinearSolver::solveWithPCG
(
    Eigen::Ref<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B
)
{
    if (x.size() != B.size())
    {
        FatalError("PCG: x and B size mismatch");
    }

    pcg_.setMaxIterations(maxIterations_);
    pcg_.setTolerance(tolerance_);
    pcg_.preconditioner().setInitialShift(icInitialShift_);

    // Sparsity analysis
    if (!solverInitialized_)
    {
        pcg_.analyzePattern(A);
        solverInitialized_ = true;
    }

    // Numeric factorization
    pcg_.factorize(A);

    if (pcg_.info() != Eigen::Success)
    {
        Warning
        (
            "Field '" + fieldName_
          + "': PCG factorization failed!"
        );
        return;
    }

    x = pcg_.solveWithGuess(B, x);

    if (debug_)
    {
        std::cout
            << "  [" << fieldName_ << "] PCG: "
            << pcg_.iterations() << " iterations"
            << ", residual = " << std::scientific
            << pcg_.error() << std::fixed
            << std::endl;
    }
}