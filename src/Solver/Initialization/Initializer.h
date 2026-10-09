/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Initializer.h
 * @brief Base class for flow field initialization
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <memory>

// Project headers
#include "StringTypes.h"
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "CellData.h"

// *************************** Forward Declarations ***************************

struct CaseConfiguration;
class GradientScheme;
class LinearSolver;

// ***************************** class Initializer ****************************

class Initializer
{
public:

// ************************* Special Member Functions *************************

    /// Default constructor
    Initializer() = default;

    /// Copy constructor and assignment - Not copyable
    Initializer(const Initializer&) = delete;
    Initializer& operator=(const Initializer&) = delete;

    /// Move constructor and assignment - Not movable
    Initializer(Initializer&&) = delete;
    Initializer& operator=(Initializer&&) = delete;

    /// Destructor
    virtual ~Initializer() = default;

// ***************************** Runtime Selection ****************************

    /// Construct the initializer selected by case configuration
    [[nodiscard]] static std::unique_ptr<Initializer> create
    (
        const CaseConfiguration& config,
        const GradientScheme& gradScheme,
        LinearSolver& pressureSolver
    );

    /// Names of every selectable initialization mode
    [[nodiscard]] static NameList availableTypes();

// *************************** Initialization Method **************************

    /// Apply the initial condition to velocity and pressure fields
    virtual void initialize
    (
        const Mesh& mesh,
        BoundaryConditions& bc,
        ScalarField& Ux,
        ScalarField& Uy,
        ScalarField& Uz,
        ScalarField& p
    ) const = 0;
};
