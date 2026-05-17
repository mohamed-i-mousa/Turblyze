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
#include <string_view>

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
    const std::ios::fmtflags flags_;
    const std::streamsize precision_;
};


namespace Logger
{
    /**
     * @brief Print a generic 80-char framed banner with the given title
     * @param title Title text shown on the second line of the banner
     */
    void sectionHeader(std::string_view title);

    /**
     * @brief Print the per-iteration banner
     * @param n Iteration number
     */
    void iterationHeader(int n);

    /// Print the closing rule that terminates a framed block
    void iterationFooter();

    /// Print the column header row for the linear solver configuration table
    void linearSolverConfigHeader();

    /**
     * @brief Print one row of the linear solver configuration table
     * @param equation Label for the equation column ("U", "p", "k", "omega")
     * @param solver Solver name ("BiCGSTAB", "PCG")
     * @param tolerance Linear solver tolerance
     * @param maxIters Linear solver iteration cap
     */
    void linearSolverConfigRow
    (
        std::string_view equation,
        std::string_view solver,
        Scalar tolerance,
        int maxIters
    );

    /**
     * @brief Print one indented "label  value" row with a Scalar value
     * @param label Row label (left-justified, padded to 24 chars)
     * @param value Scalar value, printed in scientific notation
     */
    void keyValue(std::string_view label, Scalar value);

    /**
     * @brief Print one indented "label  value" row with an int value
     * @param label Row label (left-justified, padded to 24 chars)
     * @param value Integer value, right-aligned to 12 chars
     */
    void keyValue(std::string_view label, int value);

    /**
     * @brief Print one indented "label  value" row with a string value
     * @param label Row label (left-justified, padded to 24 chars)
     * @param value String value (e.g. "kOmegaSST", "[-1e5, 1e5] Pa")
     */
    void keyValue(std::string_view label, std::string_view value);

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
        std::string_view equation,
        std::string_view solver,
        int iterations,
        Scalar linearSolverResidual
    );

    /**
     * @brief Print a sub-section title line with two-space indentation
     * @param title Title text (without trailing newline)
     */
    void subsection(std::string_view title);

    /**
     * @brief Print one min/max/mean statistics line for a scalar field
     * @param name Field label ("k", "omega", "nut")
     * @param minVal Minimum value over all cells
     * @param maxVal Maximum value over all cells
     * @param meanVal Cell-averaged mean value
     */
    void scalarStat
    (
        std::string_view name,
        Scalar minVal,
        Scalar maxVal,
        Scalar meanVal
    );

    /**
     * @brief Print one labelled scaled-residual line
     * @param name  Residual label ("mass", "velocity", ...)
     * @param value Scaled residual value
     */
    void scaledResidual(std::string_view name, Scalar value);
}
