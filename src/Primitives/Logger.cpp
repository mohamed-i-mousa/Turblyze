/******************************************************************************
 * @file Logger.cpp
 * @brief Implementations of the stateless solver-output formatting helpers
 *****************************************************************************/

#include "Logger.h"

#include <iomanip>
#include <iostream>
#include <string>


void Logger::sectionHeader(const std::string& title)
{
    StreamStateGuard guard(std::cout);
    std::cout
        << "========================================"
        << "========================================"
        << "\n " << title << "\n"
        << "----------------------------------------"
        << "----------------------------------------"
        << std::endl;
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
        << std::endl;
}


void Logger::residualTableHeader()
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "  "
        << std::left  << std::setw(11) << "Equation"
        << std::left  << std::setw(11) << "Solver"
        << std::right << std::setw(5)  << "Iters"
        << "    " << "Linear Solver Residual"
        << std::endl;

    std::cout
        << "  "
        << std::left  << std::setw(11) << "--------"
        << std::left  << std::setw(11) << "--------"
        << std::right << std::setw(5)  << "-----"
        << "    " << "-----------"
        << std::endl;
}


void Logger::residualRow
(
    const std::string& equation,
    const std::string& solver,
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
        << linearSolverResidual
        << std::endl;
}


void Logger::subsection(const std::string& title)
{
    StreamStateGuard guard(std::cout);
    std::cout << '\n' << "  " << title << std::endl;
}


void Logger::scalarStat
(
    const std::string& name,
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
        << "  mean=" << meanVal
        << std::endl;
}


void Logger::scaledResidual(const std::string& name, Scalar value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    "
        << std::left << std::setw(10) << name
        << std::scientific << std::setprecision(6)
        << value
        << std::endl;
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
        << std::endl;

    std::cout
        << "    "
        << std::left  << std::setw(11) << "--------"
        << std::left  << std::setw(11) << "--------"
        << std::left  << std::setw(12) << "----------"
        << std::right << std::setw(13) << "---------"
        << std::endl;
}


void Logger::linearSolverConfigRow
(
    const std::string& equation,
    const std::string& solver,
    Scalar tolerance,
    int maxIters
)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    "
        << std::left  << std::setw(11) << equation
        << std::left  << std::setw(11) << solver
        << std::scientific << std::setprecision(6) << tolerance
        << std::right << std::setw(13) << maxIters
        << std::endl;
}


void Logger::keyValue(const std::string& label, Scalar value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    "
        << std::left << std::setw(24) << label
        << "  "
        << std::scientific << std::setprecision(6) << value
        << std::endl;
}


void Logger::keyValue(const std::string& label, int value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    "
        << std::left  << std::setw(24) << label
        << "  "
        << std::right << std::setw(12) << value
        << std::endl;
}


void Logger::keyValue(const std::string& label, const std::string& value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    "
        << std::left << std::setw(24) << label
        << "  "
        << value
        << std::endl;
}
