/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CaseConfiguration.cpp
 * @brief Case file to typed runtime configuration mapping
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "CaseConfiguration.h"

// Standard library headers
#include <iostream>

// Project headers
#include "CaseReader.h"
#include "ConvectionSchemes.h"
#include "ErrorHandler.h"
#include "GradientScheme.h"
#include "LinearSolvers.h"
#include "RuntimeSelection.h"
#include "TurbulenceModel.h"
#include "kOmegaSST.h"

// ***************************** Internal Helpers *****************************

namespace
{

/// Fail unless selection is one of the available names for a selection point
void validateSelection
(
    const Name& selection,
    const NameList& available,
    const Name& context
)
{
    if (RuntimeSelection::isKnown(selection, available))
    {
        return;
    }

    FatalError
    (
        context + ": unknown selection '" + selection
      + "'. Valid options: " + RuntimeSelection::joinNames(available) + "."
    );
}


LinearSolverSettings momentumSolverDefaults()
{
    return LinearSolverSettings
    {
        .solver = "BiCGSTAB",
        .preconditioner = "Jacobi",
        .tolerance = S(1e-8),
        .maxIter = 1000
    };
}


LinearSolverSettings pressureSolverDefaults()
{
    return LinearSolverSettings
    {
        .solver = "PCG",
        .preconditioner = "Jacobi",
        .tolerance = S(1e-6),
        .maxIter = 1000
    };
}


LinearSolverSettings turbulenceSolverDefaults()
{
    return LinearSolverSettings
    {
        .solver = "BiCGSTAB",
        .preconditioner = "Jacobi",
        .tolerance = S(1e-6),
        .maxIter = 1000
    };
}


LinearSolverSettings readSolverEntry
(
    const CaseReader& solvers,
    const Name& key,
    LinearSolverSettings defaults
)
{
    if (solvers.hasSection(key))
    {
        const auto& entry = solvers.section(key);

        defaults.solver =
            entry.lookupOrDefault<Name>("solver", defaults.solver);
        defaults.preconditioner =
            entry.lookupOrDefault<Name>
            (
                "preconditioner",
                defaults.preconditioner
            );
        defaults.tolerance =
            entry.lookupOrDefault<Scalar>("tolerance", defaults.tolerance);
        defaults.maxIter =
            entry.lookupOrDefault<Count>("maxIter", defaults.maxIter);
    }

    validateSelection
    (
        defaults.solver,
        LinearSolver::availableSolvers(),
        "linearSolvers." + key + ".solver"
    );

    if (defaults.tolerance <= S(0.0))
    {
        FatalError("linearSolvers." + key + ".tolerance must be positive.");
    }
    if (defaults.maxIter == 0)
    {
        FatalError
        (
            "linearSolvers." + key + ".maxIter must be a positive integer."
        );
    }

    return defaults;
}


void readLinearSolvers
(
    const CaseReader& reader,
    CaseConfiguration& config
)
{
    config.linearSolvers.momentum = momentumSolverDefaults();
    config.linearSolvers.pressure = pressureSolverDefaults();
    config.linearSolvers.k = turbulenceSolverDefaults();
    config.linearSolvers.omega = turbulenceSolverDefaults();

    if (!reader.hasSection("linearSolvers"))
    {
        return;
    }

    const auto& solvers = reader.section("linearSolvers");

    config.linearSolvers.momentum =
        readSolverEntry(solvers, "U", config.linearSolvers.momentum);
    config.linearSolvers.pressure =
        readSolverEntry(solvers, "p", config.linearSolvers.pressure);

    if (!TurbulenceModel::isLaminar(config.turbulenceModel))
    {
        config.linearSolvers.k =
            readSolverEntry(solvers, "k", config.linearSolvers.k);
        config.linearSolvers.omega =
            readSolverEntry(solvers, "omega", config.linearSolvers.omega);
    }
}


void readGradientScheme
(
    const CaseReader& reader,
    CaseConfiguration& config
)
{
    const auto& schemesDict = reader.section("numericalSchemes");

    config.schemes.gradientScheme =
        schemesDict.lookupOrDefault<Name>
        (
            "gradient",
            "leastSquares"
        );

    validateSelection
    (
        config.schemes.gradientScheme,
        GradientScheme::availableSchemes(),
        "numericalSchemes.gradient"
    );
}


void readConvectionSchemes
(
    const CaseReader& reader,
    CaseConfiguration& config
)
{
    const auto& schemesDict = reader.section("numericalSchemes");

    if (!schemesDict.hasSection("convection"))
    {
        FatalError
        (
            "Missing 'convection' sub-section "
            "in numericalSchemes"
        );
    }

    const auto& convection = schemesDict.section("convection");

    config.schemes.defaultScheme =
        convection.lookupOrDefault<Name>("default", "Upwind");
    config.schemes.momentumScheme =
        convection.lookupOrDefault<Name>("U", "");
    config.schemes.kScheme =
        convection.lookupOrDefault<Name>("k", "");
    config.schemes.omegaScheme =
        convection.lookupOrDefault<Name>("omega", "");

    const NameList convectionNames = ConvectionSchemes::availableSchemes();

    validateSelection
    (
        config.schemes.defaultScheme,
        convectionNames,
        "numericalSchemes.convection.default"
    );

    if (!config.schemes.momentumScheme.empty())
    {
        validateSelection
        (
            config.schemes.momentumScheme,
            convectionNames,
            "numericalSchemes.convection.U"
        );
    }
    if (!config.schemes.kScheme.empty())
    {
        validateSelection
        (
            config.schemes.kScheme,
            convectionNames,
            "numericalSchemes.convection.k"
        );
    }
    if (!config.schemes.omegaScheme.empty())
    {
        validateSelection
        (
            config.schemes.omegaScheme,
            convectionNames,
            "numericalSchemes.convection.omega"
        );
    }

    if (config.debug)
    {
        std::cout
            << "Default convection scheme: "
            << config.schemes.defaultScheme << '\n';

        if (!config.schemes.momentumScheme.empty())
        {
            std::cout
                << "Momentum convection scheme: "
                << config.schemes.momentumScheme << '\n';
        }
        if (!config.schemes.kScheme.empty())
        {
            std::cout
                << "k convection scheme: "
                << config.schemes.kScheme << '\n';
        }
        if (!config.schemes.omegaScheme.empty())
        {
            std::cout
                << "omega convection scheme: "
                << config.schemes.omegaScheme << '\n';
        }
    }
}

} // namespace

// *************************** namespace CaseConfig ***************************

namespace CaseConfig
{

CaseConfiguration loadConfiguration(const CaseReader& reader)
{
    CaseConfiguration config;

    const auto& mesh = reader.section("mesh");
    config.meshFile = mesh.lookup<FilePath>("file");
    config.checkQuality = mesh.lookupOrDefault<bool>("checkQuality", true);

    config.numThreads = 1;
    if (reader.hasSection("parallelism"))
    {
        const auto& parallelism = reader.section("parallelism");
        const int n = parallelism.lookupOrDefault<int>("numThreads", 1);

        if (n <= 0)
        {
            Warning
            (
                "parallelism.numThreads must be a positive integer; "
                "defaulting to 1 (serial)."
            );
        }
        else
        {
            config.numThreads = static_cast<Count>(n);
        }
    }

    const auto& physicalProperties = reader.section("physicalProperties");

    config.rho = physicalProperties.lookup<Scalar>("rho");
    config.mu = physicalProperties.lookup<Scalar>("mu");
    if (config.rho <= S(0.0))
    {
        FatalError("physicalProperties.rho must be positive.");
    }
    if (config.mu <= S(0.0))
    {
        FatalError("physicalProperties.mu must be positive.");
    }

    const auto& initialConditions = reader.section("initialConditions");
    config.initialVelocity = initialConditions.lookup<Vector>("U");
    config.initialPressure = initialConditions.lookup<Scalar>("p");

    const auto& outputDict = reader.section("output");
    config.debug = outputDict.lookupOrDefault<bool>("debug", false);

    readGradientScheme(reader, config);
    readConvectionSchemes(reader, config);

    const auto& simple = reader.section("SIMPLE");

    config.maxIterations = simple.lookup<Count>("numIterations");
    if (config.maxIterations == 0)
    {
        FatalError("SIMPLE.numIterations must be a positive integer.");
    }

    config.convergenceTolerance =
        simple.lookup<Scalar>("convergenceTolerance");
    if (config.convergenceTolerance <= S(0.0))
    {
        FatalError("SIMPLE.convergenceTolerance must be positive.");
    }

    config.nNonOrthogonalCorrectors =
        simple.lookupOrDefault<Count>("nNonOrthogonalCorrectors", 0);

    const auto& relaxFactors = simple.section("relaxationFactors");
    config.alphaU = relaxFactors.lookup<Scalar>("U");
    config.alphaP = relaxFactors.lookup<Scalar>("p");
    config.alphaK = relaxFactors.lookupOrDefault<Scalar>("k", S(0.5));
    config.alphaOmega =
        relaxFactors.lookupOrDefault<Scalar>("omega", S(0.5));

    if (config.alphaU <= S(0.0) || config.alphaU > S(2.0))
    {
        FatalError("SIMPLE.relaxationFactors.U must be in (0, 2].");
    }
    if (config.alphaP <= S(0.0) || config.alphaP > S(2.0))
    {
        FatalError("SIMPLE.relaxationFactors.p must be in (0, 2].");
    }
    if (config.alphaK <= S(0.0) || config.alphaK > S(2.0))
    {
        FatalError("SIMPLE.relaxationFactors.k must be in (0, 2].");
    }
    if (config.alphaOmega <= S(0.0) || config.alphaOmega > S(2.0))
    {
        FatalError("SIMPLE.relaxationFactors.omega must be in (0, 2].");
    }

    const auto& turbulence = reader.section("turbulence");
    config.turbulenceModel = turbulence.lookup<Name>("model");

    validateSelection
    (
        config.turbulenceModel,
        TurbulenceModel::availableModels(),
        "turbulence.model"
    );

    config.turbulenceIntensity =
        turbulence.lookupOrDefault<Scalar>
        (
            "turbulenceIntensity",
            S(0.05)
        );
    config.hydraulicDiameter =
        turbulence.lookupOrDefault<Scalar>
        (
            "hydraulicDiameter",
            S(0.01)
        );

    if (config.turbulenceIntensity <= S(0.0))
    {
        FatalError("turbulence.turbulenceIntensity must be positive.");
    }
    if (config.hydraulicDiameter <= S(0.0))
    {
        FatalError("turbulence.hydraulicDiameter must be positive.");
    }

    config.initialK = S(0.0);
    config.initialOmega = S(0.0);

    if (!TurbulenceModel::isLaminar(config.turbulenceModel))
    {
        const Scalar defaultK =
            kOmegaSST::inletK
            (
                config.initialVelocity,
                config.turbulenceIntensity
            );
        const Scalar defaultOmega =
            kOmegaSST::inletOmega(defaultK, config.hydraulicDiameter);

        config.initialK =
            initialConditions.lookupOrDefault<Scalar>("k", defaultK);
        config.initialOmega =
            initialConditions.lookupOrDefault<Scalar>
            (
                "omega",
                defaultOmega
            );

        if (config.initialK < S(0.0))
        {
            FatalError("initialConditions.k must be non-negative.");
        }
        if (config.initialOmega <= S(0.0))
        {
            FatalError("initialConditions.omega must be positive.");
        }
    }

    readLinearSolvers(reader, config);

    config.vtkOutputFilename = outputDict.lookup<FilePath>("filename");
    if (config.vtkOutputFilename.empty())
    {
        FatalError("output.filename must not be empty.");
    }

    std::cout
        << "Case file loaded." << '\n';

    config.velocityConstraintEnabled = false;
    config.pressureConstraintEnabled = false;
    config.maxVelocityMagnitude = S(0.0);
    config.minPressure = S(0.0);
    config.maxPressure = S(0.0);

    if (reader.hasSection("constraints"))
    {
        const auto& constraintsDict = reader.section("constraints");

        if (constraintsDict.hasSection("velocity"))
        {
            const auto& velConstraint =
                constraintsDict.section("velocity");

            config.velocityConstraintEnabled =
                velConstraint.lookup<bool>("enabled");
            config.maxVelocityMagnitude =
                velConstraint.lookup<Scalar>("maxVelocity");

            if
            (
                config.velocityConstraintEnabled
             && config.maxVelocityMagnitude <= S(0.0)
            )
            {
                FatalError
                (
                    "constraints.velocity.maxVelocity must be positive."
                );
            }
        }

        if (constraintsDict.hasSection("pressure"))
        {
            const auto& presConstraint =
                constraintsDict.section("pressure");

            config.pressureConstraintEnabled =
                presConstraint.lookup<bool>("enabled");
            config.minPressure =
                presConstraint.lookup<Scalar>("minPressure");
            config.maxPressure =
                presConstraint.lookup<Scalar>("maxPressure");

            if
            (
                config.pressureConstraintEnabled
             && config.minPressure >= config.maxPressure
            )
            {
                FatalError
                (
                    "constraints.pressure.minPressure must be "
                    "less than maxPressure."
                );
            }
        }
    }

    config.forcesEnabled = false;
    config.forcesPatch = "";
    config.dragDirection = Vector(S(0.0), S(0.0), S(0.0));
    config.liftDirection = Vector(S(0.0), S(0.0), S(0.0));
    config.referenceVelocity = Vector(S(0.0), S(0.0), S(0.0));
    config.referenceArea = S(0.0);

    if (reader.hasSection("forces"))
    {
        const auto& forcesDict = reader.section("forces");

        config.forcesEnabled = forcesDict.lookup<bool>("enabled");

        if (config.forcesEnabled)
        {
            config.forcesPatch = forcesDict.lookup<Name>("patch");
            config.dragDirection = forcesDict.lookup<Vector>("dragDirection");
            config.liftDirection = forcesDict.lookup<Vector>("liftDirection");
            config.referenceVelocity =
                forcesDict.lookup<Vector>("referenceVelocity");
            config.referenceArea = forcesDict.lookup<Scalar>("referenceArea");

            if (config.forcesPatch.empty())
            {
                FatalError("forces.patch must not be empty.");
            }

            if (config.referenceArea <= vSmallValue)
            {
                FatalError("forces.referenceArea must be positive.");
            }

            const Scalar dragMagnitude = magnitude(config.dragDirection);
            const Scalar liftMagnitude = magnitude(config.liftDirection);
            const Scalar referenceSpeed = magnitude(config.referenceVelocity);

            if (dragMagnitude <= vSmallValue)
            {
                FatalError("forces.dragDirection must be a non-zero vector.");
            }
            if (liftMagnitude <= vSmallValue)
            {
                FatalError("forces.liftDirection must be a non-zero vector.");
            }
            if (referenceSpeed <= vSmallValue)
            {
                FatalError
                (
                    "forces.referenceVelocity must be a non-zero vector."
                );
            }

            // Store unit directions so the force module projects cleanly
            config.dragDirection = normalized(config.dragDirection);
            config.liftDirection = normalized(config.liftDirection);
        }
    }

    return config;
}

} // namespace CaseConfig
