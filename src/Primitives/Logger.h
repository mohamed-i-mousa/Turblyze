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
    /// Print a generic 80-char framed banner with the given title
    void sectionHeader(std::string_view title);

    /// Print the per-iteration banner
    void iterationHeader(int n);

    /// Print the closing rule that terminates a framed block
    void iterationFooter();

    /// Print the column header row for the linear solver configuration table
    void linearSolverConfigHeader();

    /// Print one row of the linear solver configuration table
    void linearSolverConfigRow
    (
        std::string_view equation,
        std::string_view solver,
        Scalar tolerance,
        int maxIters
    );

    /// Print one indented label-value row with a Scalar value
    void keyValue(std::string_view label, Scalar value);

    /// Print one indented label-value row with an int value
    void keyValue(std::string_view label, int value);

    /// Print one indented label-value row with a string value
    void keyValue(std::string_view label, std::string_view value);

    /// Print the column header row for the table
    void residualTableHeader();

    /// Print one row of the per-iteration residual table
    void residualRow
    (
        std::string_view equation,
        std::string_view solver,
        int iterations,
        Scalar linearSolverResidual
    );

    /// Print a sub-section title line with two-space indentation
    void subsection(std::string_view title);

    /// Print one min/max/mean statistics line for a scalar field
    void scalarStat
    (
        std::string_view name,
        Scalar minVal,
        Scalar maxVal,
        Scalar meanVal
    );

    /// Print one labelled scaled-residual line
    void scaledResidual(std::string_view name, Scalar value);

    /// Print the non-debug one-line per-iteration residual summary
    void residualSummary(Scalar mass, Scalar velocity, Scalar pressure);

    /// Print the non-debug one-line residual summary with turbulence
    void residualSummary
    (
        Scalar mass,
        Scalar velocity,
        Scalar pressure,
        Scalar k,
        Scalar omega
    );
}
