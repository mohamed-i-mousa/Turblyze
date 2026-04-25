/******************************************************************************
 * @file Logger.h
 * @brief Stateless formatting helpers for solver console output
 *
 * @details This header defines StreamStateGuard and the Logger namespace.
 * Logger provides column-aligned helpers for the solver iteration loop and
 * the k-omega SST model when debug-mode output is enabled. All helpers save
 * and restore std::cout's format flags and precision via StreamStateGuard so
 * that scientific/precision changes never leak into unrelated output.
 *
 * @class StreamStateGuard
 * - Saves std::ostream format flags and precision on construction
 * - Restores them on destruction
 * - Non-copyable and non-movable (holds an ostream reference)
 *****************************************************************************/

#pragma once

#include <ios>
#include <ostream>
#include <string>

#include "Scalar.h"


class StreamStateGuard
{
public:

    explicit StreamStateGuard(std::ostream& os) noexcept
    :
        os_{os},
        flags_{os.flags()},
        precision_{os.precision()}
    {}

    /// Copy constructor and assignment - Not copyable (T& member)
    StreamStateGuard(const StreamStateGuard&) = delete;
    StreamStateGuard& operator=(const StreamStateGuard&) = delete;

    /// Move constructor and assignment - Not movable (T& member)
    StreamStateGuard(StreamStateGuard&&) = delete;
    StreamStateGuard& operator=(StreamStateGuard&&) = delete;

    /// Destructor restores the original format flags and precision
    ~StreamStateGuard() noexcept
    {
        os_.flags(flags_);
        os_.precision(precision_);
    }

private:

    std::ostream& os_;
    std::ios::fmtflags flags_;
    std::streamsize precision_;
};


namespace Logger
{
    /**
     * @brief Print the per-iteration banner
     * @param n Iteration number
     */
    void iterationHeader(int n);

    /// Print the closing rule that terminates an iteration block
    void iterationFooter();

    /// Print the column header row for the table
    void residualTableHeader();

    /**
     * @brief Print one row of the per-iteration residual table
     * @param equation Label for the equation column ("Ux", "p'", "k", ...)
     * @param solver Solver name ("BiCGSTAB", "PCG")
     * @param iterations Linear-solver iterations performed
     * @param linearSolverResidual Final linear-solver residual
     */
    void residualRow
    (
        const std::string& equation,
        const std::string& solver,
        int iterations,
        Scalar linearSolverResidual
    );

    /**
     * @brief Print a sub-section title line with two-space indentation
     * @param title Title text (without trailing newline)
     */
    void subsection(const std::string& title);

    /**
     * @brief Print one min/max/mean statistics line for a scalar field
     * @param name Field label ("k", "omega", "nut")
     * @param minVal Minimum value over all cells
     * @param maxVal Maximum value over all cells
     * @param meanVal Cell-averaged mean value
     */
    void scalarStat
    (
        const std::string& name,
        Scalar minVal,
        Scalar maxVal,
        Scalar meanVal
    );

    /**
     * @brief Print one labelled scaled-residual line
     * @param name  Residual label ("mass", "velocity", ...)
     * @param value Scaled residual value
     */
    void scaledResidual(const std::string& name, Scalar value);
}
