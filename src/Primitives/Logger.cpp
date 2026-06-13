/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Logger.cpp
 * @brief Implementations of the stateless solver-output formatting helpers
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Logger.h"

// Standard library headers
#include <iomanip>
#include <iostream>
#include <string>

// ***************************** namespace Logger *****************************

void Logger::sectionHeader(MessageRef title)
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


void Logger::iterationHeader(Count n)
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
    NameRef equation,
    NameRef solver,
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


void Logger::subsection(MessageRef title)
{
    StreamStateGuard guard(std::cout);
    std::cout << '\n' << "  " << title << '\n';
}


void Logger::breakdownHeader(MessageRef cornerLabel)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << '\n' << "  "
        << std::left  << std::setw(14) << cornerLabel
        << std::right
        << std::setw(16) << "Pressure"
        << std::setw(16) << "Friction"
        << std::setw(16) << "Total" << '\n';

    std::cout
        << "  "
        << std::left  << std::setw(14) << "----------"
        << std::right
        << std::setw(16) << "--------"
        << std::setw(16) << "--------"
        << std::setw(16) << "-----" << '\n';
}


void Logger::breakdownRow
(
    MessageRef label,
    Scalar pressure,
    Scalar friction,
    Scalar total
)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "  "
        << std::left  << std::setw(14) << label
        << std::scientific << std::setprecision(6)
        << std::right
        << std::setw(16) << pressure
        << std::setw(16) << friction
        << std::setw(16) << total << '\n';
}


void Logger::scalarStat
(
    NameRef name,
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


void Logger::scaledResidual(NameRef name, Scalar value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    "
        << std::left << std::setw(10) << name
        << std::scientific << std::setprecision(6)
        << value << '\n';
}


void Logger::residualSummary
(
    Scalar mass,
    Scalar velocity,
    Scalar pressure
)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << " - Mass: " << std::scientific << mass
        << ", Velocity: " << velocity
        << ", Pressure: " << pressure << '\n';
}


void Logger::residualSummary
(
    Scalar mass,
    Scalar velocity,
    Scalar pressure,
    std::span<const Residuals> residuals
)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << " - Mass: " << std::scientific << mass
        << ", Velocity: " << velocity
        << ", Pressure: " << pressure;

    for (const Residuals& residual : residuals)
    {
        std::cout
            << ", " << residual.first << ": " << residual.second;
    }

    std::cout << '\n';
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
    NameRef equation,
    NameRef solver,
    Scalar tolerance,
    Count maxIters
)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left  << std::setw(11) << equation
        << std::left  << std::setw(11) << solver
        << std::scientific << std::setprecision(6) << tolerance
        << std::right << std::setw(13) << maxIters << '\n';
}


void Logger::keyValue(MessageRef label, Scalar value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left << std::setw(24) << label
        << "  " << std::scientific << std::setprecision(6) << value << '\n';
}


void Logger::keyValue(MessageRef label, Scalar value, MessageRef unit)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left << std::setw(24) << label
        << "  " << std::fixed << std::setprecision(6) << value
        << ' ' << unit << '\n';
}


void Logger::keyValue(MessageRef label, int value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left  << std::setw(24) << label
        << "  " << std::right << std::setw(12) << value << '\n';
}


void Logger::keyValue(MessageRef label, Count value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left  << std::setw(24) << label
        << "  " << std::right << std::setw(12) << value << '\n';
}


void Logger::keyValue(MessageRef label, MessageRef value)
{
    StreamStateGuard guard(std::cout);

    std::cout
        << "    " << std::left << std::setw(24) << label
        << "  " << value << '\n';
}
