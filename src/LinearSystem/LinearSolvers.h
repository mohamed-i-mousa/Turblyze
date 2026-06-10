/******************************************************************************
 * @file LinearSolvers.h
 * @brief Iterative solver hierarchy for sparse linear systems
 *
 * @details Polymorphic LinearSolver interface with an EigenLinearSolver
 * adapter that wraps an Eigen iterative solver. The Eigen solver type is
 * fixed by the EigenLinearSolver template argument; concrete classes
 * (BiCGSTAB, PCG) inherit a specific instantiation and supply the
 * algorithm name used in diagnostics.
 *
 * The shared solve path lives in EigenLinearSolver<T>::solve() and is
 * defined inline, so this header pulls Eigen's iterative-solver template
 * into every translation unit that includes it.
 *
 * @struct SolvePerformance
 * - Convergence diagnostics cached by LinearSolver after solve()
 *
 * @class LinearSolver
 * - Abstract base with virtual solve()
 * - Holds shared convergence configuration and cached solve diagnostics
 *
 * @class EigenLinearSolver
 * - Templated adapter that owns the concrete Eigen iterative solver
 * - Implements the common analyse/factorize/solveWithGuess workflow once
 * - Stores sparsity-pattern analysis state per solver instance
 *
 * @class BiCGSTAB
 * - Final solver class inheriting EigenLinearSolver<EigenBiCGSTAB> (Jacobi
 *   preconditioner); use for non-symmetric momentum/turbulence systems
 *
 * @class PCG
 * - Final solver class inheriting EigenLinearSolver<EigenPCG> (Jacobi
 *   preconditioner); use for symmetric positive definite systems (p')
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <limits>

// External library headers
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/IterativeLinearSolvers>

// Project headers
#include "ErrorHandler.h"
#include "Scalar.h"
#include "StringTypes.h"

// ********************************** Aliases *********************************

// Eigen type reductions for readability
using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;
using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
using VecRef = Eigen::Ref<Vec>;
using EigenVectorMap = Eigen::Map<Vec>;
using JacobiPreconditioner = Eigen::DiagonalPreconditioner<Scalar>;
static constexpr int LowerUpper = Eigen::Lower | Eigen::Upper;

using EigenBiCGSTAB =
    Eigen::BiCGSTAB<SparseMatrix, JacobiPreconditioner>;
using EigenPCG =
    Eigen::ConjugateGradient
    <SparseMatrix, LowerUpper, JacobiPreconditioner>;

// ************************** struct SolvePerformance *************************

struct SolvePerformance
{
    /// Solver name (view into LinearSolver::name() — program-lifetime storage)
    NameRef solverName = {};

    /// Iterations performed by the solve call
    int iterations = 0;

    /// Final relative residual reported by Eigen
    Scalar finalResidual = S(0.0);

    /// Whether Eigen reported successful convergence
    bool converged = false;
};

// **************************** class LinearSolver ****************************

class LinearSolver
{
public:

// ************************* Special Member Functions *************************

    /// Construct linear solver with convergence parameters
    LinearSolver
    (
        Scalar tolerance = S(1e-6),
        Count maxIterations = 1000
    )
    :
        tolerance_{tolerance},
        maxIterations_{maxIterations}
    {}


    /// Copy constructor and assignment - Not copyable (Polymorphic)
    LinearSolver(const LinearSolver&) = delete;
    LinearSolver& operator=(const LinearSolver&) = delete;

    /// Move constructor and assignment - Not movable
    LinearSolver(LinearSolver&&) = delete;
    LinearSolver& operator=(LinearSolver&&) = delete;

    /// Virtual destructor for polymorphic deletion
    virtual ~LinearSolver() noexcept = default;

// ****************************** Setter Methods ******************************

    /// Set relative residual tolerance
    void setTolerance(Scalar tol) noexcept
    {
        tolerance_ = tol;
    }

    /// Set maximum solver iterations
    void setMaxIterations(Count maxIter) noexcept
    {
        maxIterations_ = maxIter;
    }

// ***************************** Accessor Methods *****************************

    /// Get relative tolerance
    [[nodiscard]] Scalar tolerance() const noexcept
    {
        return tolerance_;
    }

    /// Get maximum iterations
    [[nodiscard]] Count maxIterations() const noexcept
    {
        return maxIterations_;
    }

    /// Iterations performed by the last solve call (-1 = no solve yet)
    [[nodiscard]] int lastIterations() const noexcept
    {
        return lastPerformance_.iterations;
    }

    /// Final relative residual reported by the last solve call
    [[nodiscard]] Scalar lastResidual() const noexcept
    {
        return lastPerformance_.finalResidual;
    }

    /// Full diagnostics from the last solve call
    [[nodiscard]] const SolvePerformance& lastPerformance() const noexcept
    {
        return lastPerformance_;
    }

// ******************************* Solver Method ******************************

    /// Solve sparse system using the derived algorithm
    virtual void solve
    (
        VecRef x,
        const SparseMatrix& A,
        const Vec& B
    ) = 0;

    /// Algorithm label used in diagnostic output
    [[nodiscard]] virtual NameRef name() const noexcept = 0;

// ***************************** Protected Methods ****************************

protected:

    /// Store diagnostics from the most recent solve call
    void setLastPerformance(const SolvePerformance& performance) noexcept
    {
        lastPerformance_ = performance;
    }

// ****************************** Private Members *****************************

private:

    /// Relative residual tolerance for convergence
    Scalar tolerance_;

    /// Maximum solver iterations before failure
    Count maxIterations_;

    /// Diagnostics from the most recent solve call
    SolvePerformance lastPerformance_
    {
        .solverName    = {},
        .iterations    = -1,
        .finalResidual = std::numeric_limits<Scalar>::quiet_NaN(),
        .converged     = false
    };
};


// ************************** class EigenLinearSolver *************************

template<typename EigenSolverT>
class EigenLinearSolver : public LinearSolver
{
public:

    using LinearSolver::LinearSolver;

    /// Solve sparse system using the owned Eigen iterative solver
    void solve
    (
        VecRef x,
        const SparseMatrix& A,
        const Vec& B
    ) final
    {
        if (x.size() != B.size())
        {
            FatalError("LinearSolver: x and B size mismatch");
        }

        solver_.setMaxIterations(static_cast<int>(maxIterations()));
        solver_.setTolerance(tolerance());

        if (!patternAnalyzed_)
        {
            solver_.analyzePattern(A);
            patternAnalyzed_ = true;
        }

        solver_.factorize(A);

        if (solver_.info() != Eigen::Success)
        {
            Warning(Name(name()) + ": factorization failed");

            const SolvePerformance performance
            {
                .solverName    = name(),
                .iterations    = 0,
                .finalResidual = std::numeric_limits<Scalar>::quiet_NaN(),
                .converged     = false
            };

            setLastPerformance(performance);
            return;
        }

        const Vec xPrev = x;

        x = solver_.solveWithGuess(B, x);

        if (!x.allFinite())
        {
            Warning
            (
                Name(name())
              + ": non-finite solution — rolling back to previous iterate"
            );
            x = xPrev;
        }

        const bool converged = solver_.info() == Eigen::Success;

        const SolvePerformance performance
        {
            .solverName    = name(),
            .iterations    = static_cast<int>(solver_.iterations()),
            .finalResidual = solver_.error(),
            .converged     = converged
        };

        setLastPerformance(performance);
    }

private:

    /// Owned Eigen iterative solver
    EigenSolverT solver_;

    /// Whether this instance has analysed the sparsity pattern
    bool patternAnalyzed_ = false;
};

// ****************************** class BiCGSTAB ******************************

class BiCGSTAB final : public EigenLinearSolver<EigenBiCGSTAB>
{
public:

    using EigenLinearSolver<EigenBiCGSTAB>::EigenLinearSolver;

    /// Runtime selection name
    static constexpr NameRef typeName = "BiCGSTAB";

    [[nodiscard]] NameRef name() const noexcept override
    {
        return typeName;
    }
};

// ********************************* class PCG ********************************

class PCG final : public EigenLinearSolver<EigenPCG>
{
public:

    using EigenLinearSolver<EigenPCG>::EigenLinearSolver;

    /// Runtime selection name
    static constexpr NameRef typeName = "PCG";

    [[nodiscard]] NameRef name() const noexcept override
    {
        return typeName;
    }
};
