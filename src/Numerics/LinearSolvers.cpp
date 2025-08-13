#include "LinearSolvers.h"

namespace LinearSolvers {

bool BiCGSTAB(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
              const Eigen::SparseMatrix<Scalar>& A,
              const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B,
              Scalar tolerance,
              int max_iterations,
              const std::string& fieldName) 
{
    // Resize and initialize vector x     
    if (x.size() != B.size()) {
        x.resize(B.size()); 
        x.setZero();
    }

    size_t system_size = B.size();

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> r(system_size);
    Scalar r_norm_initial_avg, r_norm_final_avg;    
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> r_initial = B - A * x;
    Scalar r_sum_initial = S(0.0);
    for (size_t i = 0; i < system_size; ++i) {
        r_sum_initial += std::abs(r_initial(i));
    }
    r_norm_initial_avg = (system_size > 0) ? (r_sum_initial / static_cast<Scalar>(system_size)) : S(0.0);
    // --- Main Solve (with ILU preconditioning) ---
    Eigen::BiCGSTAB<
        Eigen::SparseMatrix<Scalar>,
        Eigen::IncompleteLUT<Scalar>
    > bicgstab;
    bicgstab.setMaxIterations(max_iterations);
    bicgstab.setTolerance(tolerance);

    // Configure ILU (ILUT) preconditioner parameters
    // Note: higher fill-factor typically improves convergence but costs memory/time.
    bicgstab.preconditioner().setFillfactor(5);
    bicgstab.preconditioner().setDroptol(S(1e-4));

    bicgstab.compute(A);
    if (bicgstab.info() != Eigen::Success) {
        std::cerr << "Error for field '" << fieldName << "': BiCGSTAB compute (matrix decomposition/analysis) failed!" << std::endl;
        std::cerr << "  Eigen Info Code: " << bicgstab.info() << std::endl;
        // Common reasons: matrix is singular, not square, or has NaNs/Infs.
        return false;
    }

    // Use solveWithGuess, as x contains the initial guess (e.g., from previous iteration or zeros)
    x = bicgstab.solveWithGuess(B, x);

    std::cout << "\n--- Solver Statistics for Field: '" << fieldName << "' ---" << std::endl;
    std::cout << "  Iterations:     " << bicgstab.iterations() << std::endl;
    Scalar estimatedError = static_cast<Scalar>(bicgstab.error());
    std::cout << "  Estimated Error (solver reported): " << estimatedError << std::endl;

    // Abort on non-finite estimated error (divergence)
    if (!std::isfinite(static_cast<double>(estimatedError))) {
        throw std::runtime_error("  Divergence detected for field '" + fieldName +
                                 "': estimated error is not finite (" + std::to_string(estimatedError) + ")");
    }
    
    bool converged = (bicgstab.info() == Eigen::Success);
    if (converged) {
        std::cout << "  Convergence:    SUCCESS" << std::endl;
    } else {
        std::cout << "  Convergence:    FAILED (Code: " << bicgstab.info();
        // Provide more human-readable error for common cases
        if (bicgstab.info() == Eigen::NumericalIssue) {
            std::cout << " - Numerical Issue";
        } else if (bicgstab.info() == Eigen::NoConvergence) {
            std::cout << " - No Convergence (max iterations reached)";
        } else if (bicgstab.info() == Eigen::InvalidInput) {
            std::cout << " - Invalid Input";
        }
        std::cout << ")" << std::endl;
        std::cerr << "  Warning: Solver for field '" << fieldName << "' did not converge to the desired tolerance." << std::endl;
    }

    // --- Calculate and Print Actual Final Residuals ---
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> r_final = B - A * x; // Calculate actual residual: r = B - Ax
    Scalar exact_residual_norm_L2 = r_final.norm(); // L2 norm (Euclidean norm)

    Scalar r_sum_final = S(0.0);
    for (size_t i = 0; i < system_size; ++i) {
        r_sum_final += std::abs(r_final(i));
    }
    r_norm_final_avg = (system_size > 0) ? (r_sum_final / static_cast<Scalar>(system_size)) : S(0.0);

    std::cout << "  Avg Initial Abs Residual (from guess): " << r_norm_initial_avg << std::endl;
    std::cout << "  Avg Final Abs Residual:                " << r_norm_final_avg << std::endl;
    
    if (r_norm_initial_avg > S(1e-12)) { // Avoid division by zero for ratio
        Scalar r_norm_ratio = r_norm_final_avg / r_norm_initial_avg;
        std::cout << "  Residual Reduction (Avg Abs Ratio):    " << r_norm_ratio << std::endl;
    }
    std::cout << "  Final Residual L2 Norm:                " << exact_residual_norm_L2 << std::endl;
    if (!std::isfinite(static_cast<double>(r_norm_final_avg)) ||
        !std::isfinite(static_cast<double>(exact_residual_norm_L2))) {
        std::cerr << "  Divergence detected for field '" << fieldName
                  << "': residual norms are not finite." << std::endl;
        return false;
    }
    std::cout << "----------------------------------------" << std::endl;

    return converged;
}

} 