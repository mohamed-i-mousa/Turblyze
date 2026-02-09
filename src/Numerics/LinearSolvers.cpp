/******************************************************************************
 * @file LinearSolvers.cpp
 * @brief LinearSolver class implementation
 *****************************************************************************/

#include "LinearSolvers.hpp"

// ************************* Constructor & Setters ************************

LinearSolver::LinearSolver
(
    const std::string& fieldName,
    Scalar tolerance,
    int maxIterations,
    Scalar relTol,
    bool verbose
)
:   fieldName_(fieldName),
    tolerance_(tolerance),
    maxIterations_(maxIterations),
    relTol_(relTol),
    verbose_(verbose)
{}

// ***************************** Shared Helpers *****************************

Scalar LinearSolver::rmsResidual
(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& 
    residualVector
) const
{
    size_t n = residualVector.size();
    return std::sqrt(residualVector.squaredNorm() / S(n));
}

void LinearSolver::reportConvergence
(
    int iterations,
    Eigen::ComputationInfo eigenInfo,
    const std::optional<ResidualInfo>& residuals
) const
{
    std::cout
        << "\n--- Field '" << fieldName_
        << "' ---" << std::endl;

    std::cout
        << "  Iterations:       " << iterations
        << std::endl;

    if (eigenInfo == Eigen::Success)
    {
        std::cout
            << "  Eigen status:     converged"
            << std::endl;
    }
    else if (eigenInfo == Eigen::NoConvergence)
    {
        std::cout
            << "  Eigen status:     max iterations"
            << std::endl;
    }
    else
    {
        std::cout
            << "  Eigen status:     failed (code "
            << eigenInfo << ")" << std::endl;
    }

    if (!residuals.has_value())
    {
        std::cout
            << "  (Exact residuals: disabled)"
            << std::endl;
        return;
    }

    std::cout
        << std::scientific
        << "  Initial residual: "
        << residuals->initialResidual << std::endl;

    std::cout
        << "  Final residual:   "
        << residuals->finalResidual;

    if (residuals->finalResidual < tolerance_)
    {
        std::cout
            << "  CONVERGED";
    }

    std::cout
        << std::endl;

    std::cout
        << "  Tolerance:        "
        << tolerance_ << std::endl;

    if
    (
        relTol_ > S(0.0)
     && residuals->initialResidual > S(1e-12)
    )
    {
        std::cout
            << "  Residual ratio:   "
            << residuals->residualRatio
            << "  (relTol: " << relTol_ << ")"
            << std::endl;
    }
}

// ***************************** BiCGSTAB Solver *****************************

bool LinearSolver::solveWithBiCGSTAB
(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B
)
{
    // Resize and initialize if needed
    if (x.size() != B.size())
    {
        x.resize(B.size());
        x.setZero();
    }

    // Compute initial residual
    Scalar initialResidual = S(0.0);
    if (computeResiduals_)
    {
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
        initialResidualVector = B - A * x;

        initialResidual =
            rmsResidual(initialResidualVector);
    }

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

    // Compute residuals if requested
    std::optional<ResidualInfo> residuals;

    if (computeResiduals_)
    {
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
        finalResVec = B - A * x;

        ResidualInfo info;

        info.initialResidual = initialResidual;

        info.finalResidual = rmsResidual(finalResVec);
        info.residualRatio =
            (initialResidual > smallValue)
            ? info.finalResidual/initialResidual
            : S(0.0);

        residuals = info;
    }

    if (verbose_)
    {
        reportConvergence
        (
            bicgstab.iterations(),
            bicgstab.info(),
            residuals
        );
    }

    return true;
}

// ******************************* PCG Solver *******************************

bool LinearSolver::solveWithPCG
(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B
)
{
    // Resize and initialize if needed
    if (x.size() != B.size())
    {
        x.resize(B.size());
        x.setZero();
    }

    // Compute initial residual only if requested
    Scalar initialResidual = S(0.0);
    if (computeResiduals_)
    {
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
        initialResidualVector = B - A * x;

        initialResidual =
            rmsResidual(initialResidualVector);
    }

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

    // Compute residuals if requested
    std::optional<ResidualInfo> residuals;

    if (computeResiduals_)
    {
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
        finalResVec = B - A * x;

        ResidualInfo info;
        info.initialResidual = initialResidual;
        info.finalResidual = rmsResidual(finalResVec);

        info.residualRatio =
            (initialResidual > smallValue)
            ? info.finalResidual/initialResidual
            : S(0.0);

        residuals = info;
    }

    if (verbose_)
    {
        reportConvergence
        (
            pcg.iterations(),
            pcg.info(),
            residuals
        );
    }

    return true;
}