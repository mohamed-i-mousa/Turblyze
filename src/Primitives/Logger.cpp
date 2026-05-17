/******************************************************************************
 * @file Logger.cpp
 * @brief Implementations of the stateless solver-output formatting helpers
 *****************************************************************************/

#include "Logger.h"

#include <iomanip>
#include <iostream>
#include <string>


void Logger::sectionHeader(std::string_view title)
{
    StreamStateGuard guard(std::cout);
    std::cout
        << "========================================"
        << "========================================"
        << '\n' << " " << title << '\n'
        << "----------------------------------------"
        << "----------------------------------------"
        << '\n';
}


void Logger::iterationHeader(int n)
{
    sectionHeader("Iteration " + std::to_string(n));
}


void Logger::iterationFooter()
{
    std::cout
        << "========================================"
        << "========================================"
        << '\n';
}


void Logger::residualTableHeader()
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "  "
        << std::left  << std::setw(11) << "Equation"
        << std::left  << std::setw(11) << "Solver"
        << std::right << std::setw(5)  << "Iters"
        << "    " << "Linear Solver Residual" << '\n';

    std::cout
        << "  "
        << std::left  << std::setw(11) << "--------"
        << std::left  << std::setw(11) << "--------"
        << std::right << std::setw(5)  << "-----"
        << "    " << "-----------" << '\n';
}


void Logger::residualRow
(
    std::string_view equation,
    std::string_view solver,
    int iterations,
    Scalar linearSolverResidual
)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "  "
        << std::left  << std::setw(11) << equation
        << std::left  << std::setw(11) << solver
        << std::right << std::setw(5)  << iterations
        << "    "
        << std::scientific << std::setprecision(6)
        << linearSolverResidual << '\n';
}


void Logger::subsection(std::string_view title)
{
    StreamStateGuard guard(std::cout);
    std::cout << '\n' << "  " << title << '\n';
}


void Logger::scalarStat
(
    std::string_view name,
    Scalar minVal,
    Scalar maxVal,
    Scalar meanVal
)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    "
        << std::left << std::setw(7) << name
        << std::scientific << std::setprecision(2)
        << "min="  << minVal
        << "  max="  << maxVal
        << "  mean=" << meanVal << '\n';
}


void Logger::scaledResidual(std::string_view name, Scalar value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    "
        << std::left << std::setw(10) << name
        << std::scientific << std::setprecision(6)
        << value << '\n';
}


void Logger::linearSolverConfigHeader()
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    "
        << std::left  << std::setw(11) << "Equation"
        << std::left  << std::setw(11) << "Solver"
        << std::left  << std::setw(12) << "Tolerance"
        << std::right << std::setw(13) << "Max Iters"
        << '\n';

    std::cout
        << "    "
        << std::left  << std::setw(11) << "--------"
        << std::left  << std::setw(11) << "--------"
        << std::left  << std::setw(12) << "----------"
        << std::right << std::setw(13) << "---------"
        << '\n';
}


void Logger::linearSolverConfigRow
(
    std::string_view equation,
    std::string_view solver,
    Scalar tolerance,
    int maxIters
)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left  << std::setw(11) << equation
        << std::left  << std::setw(11) << solver
        << std::scientific << std::setprecision(6) << tolerance
        << std::right << std::setw(13) << maxIters << '\n';
}


void Logger::keyValue(std::string_view label, Scalar value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left << std::setw(24) << label
        << "  " << std::scientific << std::setprecision(6) << value << '\n';
}


void Logger::keyValue(std::string_view label, int value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left  << std::setw(24) << label
        << "  " << std::right << std::setw(12) << value << '\n';
}


void Logger::keyValue(std::string_view label, std::string_view value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left << std::setw(24) << label
        << "  " << value << '\n';
}
