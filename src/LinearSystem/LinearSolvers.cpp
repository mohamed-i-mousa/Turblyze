/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LinearSolvers.cpp
 * @brief Runtime selection of linear solvers
 *****************************************************************************/

// ********************************** Headers *********************************

#include "LinearSolvers.h"

// Project headers
#include "RuntimeSelection.h"

// **************************** Runtime Selection ****************************

std::unique_ptr<LinearSolver> LinearSolver::create
(
    NameRef solverName,
    Scalar tolerance,
    Count maxIterations
)
{
    if (solverName == "BiCGSTAB")
    {
        return std::make_unique<BiCGSTAB>(solverName, tolerance, maxIterations);
    }

    if (solverName == "PCG")
    {
        return std::make_unique<PCG>(solverName, tolerance, maxIterations);
    }

    RuntimeSelection::unknownSelection
    (
        "linear solver",
        Name(solverName),
        availableSolvers()
    );
}


NameList LinearSolver::availableSolvers()
{
    return {"BiCGSTAB", "PCG"};
}
