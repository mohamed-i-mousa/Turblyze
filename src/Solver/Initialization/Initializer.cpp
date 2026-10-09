/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Initializer.cpp
 * @brief Runtime selection and factory for flow initializers
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Initializer.h"

// Standard library headers
#include <memory>

// Project headers
#include "RuntimeSelection.h"
#include "CaseConfiguration.h"
#include "Uniform.h"
#include "PotentialFlow.h"

// ***************************** Runtime Selection ****************************

std::unique_ptr<Initializer> Initializer::create
(
    const CaseConfiguration& config,
    const GradientScheme& gradScheme,
    LinearSolver& pressureSolver
)
{
    const Name& typeName = config.initializationType;

    if (typeName == "Uniform")
    {
        return std::make_unique<Initialization::Uniform>
        (
            config.initialVelocity,
            config.initialPressure
        );
    }


    if (typeName == "potentialFlow")
    {
        return std::make_unique<Initialization::PotentialFlow>
        (
            config.initialVelocity,
            config.initialPressure,
            config.rho,
            gradScheme,
            config.potentialFlowCorrectors,
            pressureSolver
        );
    }

    RuntimeSelection::unknownSelection
    (
        "initialization type",
        typeName,
        availableTypes()
    );
}


NameList Initializer::availableTypes()
{
    return {"Uniform", "potentialFlow"};
}
