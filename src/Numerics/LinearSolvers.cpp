/******************************************************************************
 * @file LinearSolvers.cpp
 * @brief Linear solver implementations for sparse matrix systems
 *****************************************************************************/

#include "LinearSolvers.h"

namespace LinearSolvers {

bool BiCGSTAB
(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B,
    Scalar tolerance,
    int maxIterations,
    const std::string& fieldName
) 
{
    // Resize and initialize vector x     
    if (x.size() != B.size()) 
    {
        x.resize(B.size()); 
        x.setZero();
    }

    size_t systemSize = B.size();

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> r(systemSize);
    Scalar initialResidual, finalResidual;
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> initialResidualVector = B - A * x;
    Scalar initialResidualSum = S(0.0);

    for (size_t i = 0; i < systemSize; ++i) 
    {
        initialResidualSum += std::abs(initialResidualVector(i));
    }
    initialResidual = initialResidualSum / S(systemSize);

    // --- Main Solve (with ILU preconditioning) ---
    Eigen::BiCGSTAB
    <
        Eigen::SparseMatrix<Scalar>,
        Eigen::IncompleteLUT<Scalar>
    > bicgstab;

    bicgstab.setMaxIterations(maxIterations);
    bicgstab.setTolerance(tolerance);

    // Configure ILU (ILUT) preconditioner parameters
    // higher fill-factor typically improves convergence but costs memory/time.
    bicgstab.preconditioner().setFillfactor(5);
    bicgstab.preconditioner().setDroptol(S(1e-6));

    bicgstab.compute(A);
    if (bicgstab.info() != Eigen::Success) 
    {
        std::cerr   << "Error for field '" << fieldName 
                    << "': BiCGSTAB compute failed!" << std::endl;

        std::cerr   << "  Eigen Info Code: " << bicgstab.info() << std::endl;

        return false;
    }

    // Use solveWithGuess, as x contains the initial guess 
    // (e.g., from previous iteration or zeros)

    x = bicgstab.solveWithGuess(B, x);

    std::cout   << "\n--- Solver Statistics for Field: '" 
                << fieldName << "' ---" << std::endl;

    std::cout   << "  Iterations:     " << bicgstab.iterations() << std::endl;

    Scalar estimatedError = S(bicgstab.error());

    std::cout   << "  Estimated Error (solver reported): " 
                << std::scientific << estimatedError << std::endl;

    // Abort on non-finite estimated error or excessively large error (divergence)
    if (!std::isfinite(static_cast<double>(estimatedError)) || estimatedError > 1e6) 
    {
        throw std::runtime_error
        (
            "  Divergence detected for field '" 
          + fieldName 
          + "': estimated error is not finite or too large (" 
          + std::to_string(estimatedError) 
          + ")"
        );
    }
    
    bool converged = (bicgstab.info() == Eigen::Success);
    if (converged) 
    {
        std::cout << "  Convergence:    SUCCESS" << std::endl;
    }
    else
    {
        std::cout << "  Convergence:    FAILED (Code: " << bicgstab.info();
        // Provide more human-readable error for common cases
        if (bicgstab.info() == Eigen::NumericalIssue) 
        {
            std::cout << " - Numerical Issue";
        }
        else if (bicgstab.info() == Eigen::NoConvergence) 
        {
            std::cout << " - No Convergence (max iterations reached)";
        }
        else if (bicgstab.info() == Eigen::InvalidInput)
        {
            std::cout << " - Invalid Input";
        }
        
        std::cout   << ")" << std::endl;
        
        std::cerr   << "  Warning: Solver for field '" << fieldName 
                    << "' did not converge to the desired tolerance." 
                    << std::endl;
    }

    // --- Calculate and Print Actual Final Residuals ---
    
    // Calculate actual residual: r = B - Ax
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> finalResidualVector = B - A * x;

    // L2 norm (Euclidean norm)
    Scalar exactResidualNormL2 = finalResidualVector.norm(); 

    Scalar finalResidualSum = S(0.0);
    
    for (size_t i = 0; i < systemSize; ++i) 
    {
        finalResidualSum += std::abs(finalResidualVector(i));
    }

    finalResidual = finalResidualSum / S(systemSize);

    std::cout   << "  Avg Initial Abs Residual (from guess): "
                << std::scientific << initialResidual << std::endl;

    std::cout   << "  Avg Final Abs Residual:                "
                << std::scientific << finalResidual << std::endl;
    
    if (initialResidual > S(1e-12)) 
    { // Avoid division by zero for ratio
        Scalar residualRatio = finalResidual / initialResidual;
        std::cout   << "  Residual Reduction (Avg Abs Ratio):    "
                    << std::scientific << residualRatio << std::endl;
    }
    
    std::cout   << "  Final Residual L2 Norm:                " 
                << std::scientific << exactResidualNormL2 << std::endl;

    if 
    (
        !std::isfinite(static_cast<double>(finalResidual)) 
     || !std::isfinite(static_cast<double>(exactResidualNormL2))) 
    {
        std::cerr << "  Divergence detected for field '" << fieldName
                  << "': residual norms are not finite." << std::endl;
        return false;
    }

    std::cout << "----------------------------------------" << std::endl;

    return converged;
}

bool PCG
(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
    const Eigen::SparseMatrix<Scalar>& A,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B,
    Scalar tolerance,
    int maxIterations,
    const std::string& fieldName
)
{
    // Resize and initialize vector x     
    if (x.size() != B.size()) 
    {
        x.resize(B.size()); 
        x.setZero();
    }

    size_t systemSize = B.size();

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> r(systemSize);
    Scalar initialResidual, finalResidual;
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> initialResidualVector = B - A * x;
    Scalar initialResidualSum = S(0.0);

    for (size_t i = 0; i < systemSize; ++i) 
    {
        initialResidualSum += std::abs(initialResidualVector(i));
    }
    initialResidual = initialResidualSum / S(systemSize);

    // --- Main Solve (with Incomplete Cholesky preconditioning) ---
    Eigen::ConjugateGradient
    <
        Eigen::SparseMatrix<Scalar>,
        Eigen::Lower|Eigen::Upper,
        Eigen::IncompleteCholesky<Scalar>
    > pcg;

    pcg.setMaxIterations(maxIterations);
    pcg.setTolerance(tolerance);

    // Configure Incomplete Cholesky preconditioner parameters
    // Lower initial fill-in factor for symmetric matrices
    pcg.preconditioner().setInitialShift(S(1e-3));

    pcg.compute(A);
    if (pcg.info() != Eigen::Success) 
    {
        std::cerr   << "Error for field '" << fieldName 
                    << "': PCG compute failed!" << std::endl;

        std::cerr   << "  Eigen Info Code: " << pcg.info() << std::endl;

        return false;
    }

    // Use solveWithGuess, as x contains the initial guess 
    // (e.g., from previous iteration or zeros)

    x = pcg.solveWithGuess(B, x);

    std::cout   << "\n--- Solver Statistics for Field: '" 
                << fieldName << "' ---" << std::endl;

    std::cout   << "  Iterations:     " << pcg.iterations() << std::endl;

    Scalar estimatedError = S(pcg.error());

    std::cout   << "  Estimated Error (solver reported): " 
                << std::scientific << estimatedError << std::endl;

    // Abort on non-finite estimated error or excessively large error (divergence)
    if (!std::isfinite(static_cast<double>(estimatedError)) || estimatedError > 1e6) 
    {
        throw std::runtime_error
        (
            "  Divergence detected for field '" 
          + fieldName 
          + "': estimated error is not finite or too large (" 
          + std::to_string(estimatedError) 
          + ")"
        );
    }
    
    bool converged = (pcg.info() == Eigen::Success);
    if (converged) 
    {
        std::cout << "  Convergence:    SUCCESS" << std::endl;
    }
    else
    {
        std::cout << "  Convergence:    FAILED (Code: " << pcg.info();
        // Provide more human-readable error for common cases
        if (pcg.info() == Eigen::NumericalIssue) 
        {
            std::cout << " - Numerical Issue";
        }
        else if (pcg.info() == Eigen::NoConvergence) 
        {
            std::cout << " - No Convergence (max iterations reached)";
        }
        else if (pcg.info() == Eigen::InvalidInput)
        {
            std::cout << " - Invalid Input";
        }
        
        std::cout   << ")" << std::endl;
        
        std::cerr   << "  Warning: Solver for field '" << fieldName 
                    << "' did not converge to the desired tolerance." 
                    << std::endl;
    }

    // --- Calculate and Print Actual Final Residuals ---
    
    // Calculate actual residual: r = B - Ax
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> finalResidualVector = B - A * x;

    // L2 norm (Euclidean norm)
    Scalar exactResidualNormL2 = finalResidualVector.norm(); 

    Scalar finalResidualSum = S(0.0);
    
    for (size_t i = 0; i < systemSize; ++i) 
    {
        finalResidualSum += std::abs(finalResidualVector(i));
    }

    finalResidual = finalResidualSum / S(systemSize);

    std::cout   << "  Avg Initial Abs Residual (from guess): "
                << std::scientific << initialResidual << std::endl;

    std::cout   << "  Avg Final Abs Residual:                "
                << std::scientific << finalResidual << std::endl;
    
    if (initialResidual > S(1e-12)) 
    { // Avoid division by zero for ratio
        Scalar residualRatio = finalResidual / initialResidual;
        std::cout   << "  Residual Reduction (Avg Abs Ratio):    "
                    << std::scientific << residualRatio << std::endl;
    }
    
    std::cout   << "  Final Residual L2 Norm:                " 
                << std::scientific << exactResidualNormL2 << std::endl;

    if 
    (
        !std::isfinite(static_cast<double>(finalResidual)) 
     || !std::isfinite(static_cast<double>(exactResidualNormL2))) 
    {
        std::cerr << "  Divergence detected for field '" << fieldName
                  << "': residual norms are not finite." << std::endl;
        return false;
    }

    std::cout << "----------------------------------------" << std::endl;

    return converged;
}

} 