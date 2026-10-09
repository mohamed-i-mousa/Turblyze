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
#include <format>
#include <iostream>

// Project headers
#include "CaseReader.h"
#include "ConvectionScheme.h"
#include "ErrorHandler.h"
#include "GradientScheme.h"
#include "LinearSolvers.h"
#include "MomentumTransport.h"
#include "Initializer.h"
#include "RuntimeSelection.h"
#include "TimeScheme.h"
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

    validateSelection
    (
        defaults.preconditioner,
        LinearSolver::availablePreconditioners(),
        "linearSolvers." + key + ".preconditioner"
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
    config.linearSolvers.momentum = LinearSolverSettings
    {
        .solver = "BiCGSTAB",
        .preconditioner = "Jacobi",
        .tolerance = S(1e-8),
        .maxIter = 1000
    };
    config.linearSolvers.pressure = LinearSolverSettings
    {
        .solver = "PCG",
        .preconditioner = "Jacobi",
        .tolerance = S(1e-6),
        .maxIter = 1000
    };
    config.linearSolvers.k = turbulenceSolverDefaults();
    config.linearSolvers.omega = turbulenceSolverDefaults();
    config.linearSolvers.petscOptions = "";

    if (!reader.hasSection("linearSolvers"))
    {
        return;
    }

    const auto& solvers = reader.section("linearSolvers");

    config.linearSolvers.petscOptions =
        solvers.lookupOrDefault<Name>("petscOptions", "");

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

    const NameList convectionNames = ConvectionScheme::availableSchemes();

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
        std::cout << std::format
        (
            "Default convection scheme: {}\n",
            config.schemes.defaultScheme
        );

        if (!config.schemes.momentumScheme.empty())
        {
            std::cout << std::format
            (
                "Momentum convection scheme: {}\n",
                config.schemes.momentumScheme
            );
        }
        if (!config.schemes.kScheme.empty())
        {
            std::cout << std::format
            (
                "k convection scheme: {}\n",
                config.schemes.kScheme
            );
        }
        if (!config.schemes.omegaScheme.empty())
        {
            std::cout << std::format
            (
                "omega convection scheme: {}\n",
                config.schemes.omegaScheme
            );
        }
    }
}

void readTimeControl
(
    const CaseReader& reader,
    CaseConfiguration& config
)
{
    // Defaults for the steadyState case (transient keys stay unused)
    config.time.timeStep = S(0.0);
    config.time.totalTime = S(0.0);
    config.time.writingIntervals = 1;
    config.time.nOuterCorrectors = 1;
    config.time.CrankNicolsonCoeff = S(1.0);

    const auto& time = reader.section("time");

    config.time.timeScheme = time.lookup<Name>("timeScheme");

    validateSelection
    (
        config.time.timeScheme,
        TimeScheme::availableSchemes(),
        "time.timeScheme"
    );

    if (config.time.timeScheme == "steadyState")
    {
        return;
    }

    config.time.timeStep = time.lookup<Scalar>("timeStep");
    config.time.totalTime = time.lookup<Scalar>("totalTime");
    config.time.writingIntervals =
        time.lookupOrDefault<Count>("writingIntervals", 1);
    config.time.nOuterCorrectors =
        time.lookupOrDefault<Count>("nOuterCorrectors", 1);
    config.time.CrankNicolsonCoeff =
        time.lookupOrDefault<Scalar>("CrankNicolsonCoeff", S(1.0));

    if (config.time.timeStep <= S(0.0))
    {
        FatalError("time.timeStep must be positive.");
    }
    if (config.time.totalTime <= S(0.0))
    {
        FatalError("time.totalTime must be positive.");
    }
    if (config.time.writingIntervals == 0)
    {
        FatalError("time.writingIntervals must be a positive integer.");
    }
    if (config.time.nOuterCorrectors == 0)
    {
        FatalError("time.nOuterCorrectors must be a positive integer.");
    }
    if (config.time.CrankNicolsonCoeff < S(0.0)
     || config.time.CrankNicolsonCoeff > S(1.0))
    {
        FatalError("time.CrankNicolsonCoeff must be in [0, 1].");
    }
}


void readVelocityCoupling
(
    const CaseReader& reader,
    CaseConfiguration& config
)
{
    // Default if there's no velocityCoupling section
    config.algorithm = "SIMPLE";
    config.nPrimeCorrectors = 1;

    if (reader.hasSection("velocityCoupling"))
    {
        const auto& velocityCoupling = reader.section("velocityCoupling");

        config.algorithm =
            velocityCoupling.lookupOrDefault<Name>("algorithm", "SIMPLE");

        validateSelection
        (
            config.algorithm,
            MomentumTransport::availableAlgorithms(),
            "velocityCoupling.algorithm"
        );
    }

    const bool transient = config.time.timeScheme != "steadyState";

    if (transient && config.algorithm == "SIMPLE")
    {
        Warning
        (
            "a transient time scheme cannot proceed with SIMPLE; "
            "switching algorithm to PISO."
        );
        config.algorithm = "PISO";
    }
    else if (!transient && config.algorithm == "PISO")
    {
        Warning
        (
            "the steadyState time scheme requires SIMPLE; "
            "switching algorithm to SIMPLE."
        );
        config.algorithm = "SIMPLE";
    }
}

void readSimpleControls
(
    const CaseReader& reader,
    CaseConfiguration& config
)
{
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
}

void readPisoControls
(
    const CaseReader& reader,
    CaseConfiguration& config
)
{
    const auto& piso = reader.section("PISO");

    // PISO never runs the steady solve(); maxIterations is unused.
    config.maxIterations = 1;

    config.convergenceTolerance =
        piso.lookup<Scalar>("convergenceTolerance");
    if (config.convergenceTolerance <= S(0.0))
    {
        FatalError("PISO.convergenceTolerance must be positive.");
    }

    config.nNonOrthogonalCorrectors =
        piso.lookupOrDefault<Count>("nNonOrthogonalCorrectors", 0);

    config.nPrimeCorrectors =
        piso.lookupOrDefault<Count>("nPrimeCorrectors", 1);
    if (config.nPrimeCorrectors == 0)
    {
        FatalError("PISO.nPrimeCorrectors must be a positive integer.");
    }

    // PISO is stabilized by the V/dt diagonal, not relaxation factors
    config.alphaU = S(1.0);
    config.alphaP = S(1.0);
    config.alphaK = S(1.0);
    config.alphaOmega = S(1.0);

    if (piso.hasSection("relaxationFactors"))
    {
        const auto& relaxFactors = piso.section("relaxationFactors");
        config.alphaU = relaxFactors.lookupOrDefault<Scalar>("U", S(1.0));
        config.alphaP = relaxFactors.lookupOrDefault<Scalar>("p", S(1.0));
        config.alphaK = relaxFactors.lookupOrDefault<Scalar>("k", S(1.0));
        config.alphaOmega =
            relaxFactors.lookupOrDefault<Scalar>("omega", S(1.0));
    }

    if (config.alphaU <= S(0.0) || config.alphaU > S(2.0))
    {
        FatalError("PISO.relaxationFactors.U must be in (0, 2].");
    }
    if (config.alphaP <= S(0.0) || config.alphaP > S(2.0))
    {
        FatalError("PISO.relaxationFactors.p must be in (0, 2].");
    }
    if (config.alphaK <= S(0.0) || config.alphaK > S(2.0))
    {
        FatalError("PISO.relaxationFactors.k must be in (0, 2].");
    }
    if (config.alphaOmega <= S(0.0) || config.alphaOmega > S(2.0))
    {
        FatalError("PISO.relaxationFactors.omega must be in (0, 2].");
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
    config.initializationType =
        initialConditions.lookupOrDefault<Name>("type", "Uniform");

    validateSelection
    (
        config.initializationType,
        Initializer::availableTypes(),
        "initialConditions.type"
    );

    if (config.initializationType == "Uniform")
    {
        config.initialVelocity = initialConditions.lookup<Vector>("U");
        config.initialPressure = initialConditions.lookup<Scalar>("p");
    }
    else if (config.initializationType == "potentialFlow")
    {
        config.initialVelocity = initialConditions.lookup<Vector>("U");
        config.initialPressure = initialConditions.lookup<Scalar>("p");
    }

    const auto& outputDict = reader.section("output");
    config.debug = outputDict.lookupOrDefault<bool>("debug", false);

    readGradientScheme(reader, config);
    readConvectionSchemes(reader, config);

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

    config.roughWall =
        turbulence.lookupOrDefault<bool>("roughWall", false);

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

    readTimeControl(reader, config);

    readVelocityCoupling(reader, config);

    if (config.algorithm == "SIMPLE")
    {
        readSimpleControls(reader, config);
    }
    else if (config.algorithm == "PISO")
    {
        readPisoControls(reader, config);
    }

    if (config.initializationType == "potentialFlow")
    {
        // 1 base solve + N non-orthogonal correction loops
        config.potentialFlowCorrectors = 1 + config.nNonOrthogonalCorrectors;
    }

    config.vtkOutputFilename = outputDict.lookup<FilePath>("filename");
    if (config.vtkOutputFilename.empty())
    {
        FatalError("output.filename must not be empty.");
    }

    std::cout
        << "Case file loaded." << '\n';

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
