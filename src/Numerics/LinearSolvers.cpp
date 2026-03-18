/******************************************************************************
 * @file LinearSolvers.cpp
 * @brief LinearSolver class implementation
 *****************************************************************************/

#include "LinearSolvers.hpp"

#include <cassert>
#include <iostream>

#include <eigen3/Eigen/IterativeLinearSolvers>

// ************************* Constructor & Setters ************************

LinearSolver::LinearSolver
(
    std::string fieldName,
    Scalar tolerance,
    int maxIterations
)
:   fieldName_(std::move(fieldName)),
    tolerance_(tolerance),
    maxIterations_(maxIterations)
{}


// ***************************** BiCGSTAB Solver *****************************

bool LinearSolver::solveWithBiCGSTAB
(
    Eigen::Ref<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B
)
{
    assert(x.size() == B.size());

    // Configure Eigen BiCGSTAB with ILUT preconditioning
    Eigen::BiCGSTAB
    <Eigen::SparseMatrix<Scalar>, Eigen::IncompleteLUT<Scalar>>
    bicgstab;

    bicgstab.setMaxIterations(maxIterations_);
    bicgstab.setTolerance(tolerance_);
    bicgstab.preconditioner().setFillfactor(ilutFillFactor_);
    bicgstab.preconditioner().setDroptol(ilutDropTol_);

    bicgstab.compute(A);
    if (bicgstab.info() != Eigen::Success)
    {
        std::cerr
            << "Error for field '" << fieldName_
            << "': BiCGSTAB compute failed!"
            << std::endl;

        std::cerr
            << "  Eigen Info Code: "
            << bicgstab.info() << std::endl;
        return false;
    }

    x = bicgstab.solveWithGuess(B, x);

    return true;
}

// ******************************* PCG Solver *******************************

bool LinearSolver::solveWithPCG
(
    Eigen::Ref<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B
)
{
    assert(x.size() == B.size());

    // Configure Eigen CG with Incomplete Cholesky preconditioning
    Eigen::ConjugateGradient
    <Eigen::SparseMatrix<Scalar>, Eigen::Lower|Eigen::Upper,
    Eigen::IncompleteCholesky<Scalar>>
    pcg;

    pcg.setMaxIterations(maxIterations_);
    pcg.setTolerance(tolerance_);
    pcg.preconditioner().setInitialShift(icInitialShift_);

    pcg.compute(A);
    if (pcg.info() != Eigen::Success)
    {
        std::cerr
            << "Error for field '" << fieldName_
            << "': PCG compute failed!"
            << std::endl;

        std::cerr
            << "  Eigen Info Code: "
            << pcg.info() << std::endl;
        return false;
    }

    x = pcg.solveWithGuess(B, x);

    return true;
}