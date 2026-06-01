/******************************************************************************
 * @file CaseConfiguration.cpp
 * @brief Case file to typed runtime configuration mapping
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation headers
#include "CaseConfiguration.h"

// Standard library headers
#include <iostream>

// Project headers
#include "CaseReader.h"
#include "ErrorHandler.h"
#include "kOmegaSST.h"

// ***************************** Internal helpers *****************************

namespace
{

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
    const std::string& key,
    LinearSolverSettings defaults
)
{
    if (solvers.hasSection(key))
    {
        const auto& entry = solvers.section(key);

        defaults.solver =
            entry.lookupOrDefault<std::string>("solver", defaults.solver);
        defaults.preconditioner =
            entry.lookupOrDefault<std::string>
            (
                "preconditioner",
                defaults.preconditioner
            );
        defaults.tolerance =
            entry.lookupOrDefault<Scalar>("tolerance", defaults.tolerance);
        defaults.maxIter =
            entry.lookupOrDefault<int>("maxIter", defaults.maxIter);
    }

    if (defaults.tolerance <= S(0.0))
    {
        FatalError("linearSolvers." + key + ".tolerance must be positive.");
    }
    if (defaults.maxIter <= 0)
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

    if (config.turbulenceEnabled)
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

    config.schemes.gradientName =
        schemesDict.lookupOrDefault<std::string>
        (
            "gradient",
            "leastSquares"
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

    config.schemes.defaultName =
        convection.lookupOrDefault<std::string>("default", "Upwind");
    config.schemes.momentumName =
        convection.lookupOrDefault<std::string>("U", "");
    config.schemes.kName =
        convection.lookupOrDefault<std::string>("k", "");
    config.schemes.omegaName =
        convection.lookupOrDefault<std::string>("omega", "");

    if (config.debug)
    {
        std::cout
            << "Default convection scheme: "
            << config.schemes.defaultName << '\n';

        if (!config.schemes.momentumName.empty())
        {
            std::cout
                << "Momentum convection scheme: "
                << config.schemes.momentumName << '\n';
        }
        if (!config.schemes.kName.empty())
        {
            std::cout
                << "k convection scheme: "
                << config.schemes.kName << '\n';
        }
        if (!config.schemes.omegaName.empty())
        {
            std::cout
                << "omega convection scheme: "
                << config.schemes.omegaName << '\n';
        }
    }
}

} // namespace

// **************************** namespace CaseConfig *************************

namespace CaseConfig
{

CaseConfiguration loadConfiguration(const CaseReader& reader)
{
    CaseConfiguration config;

    const auto& mesh = reader.section("mesh");
    config.meshFilePath = mesh.lookup<std::string>("file");
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
            config.numThreads = n;
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

    config.maxIterations = simple.lookup<int>("numIterations");
    if (config.maxIterations <= 0)
    {
        FatalError("SIMPLE.numIterations must be a positive integer.");
    }

    config.convergenceTolerance =
        simple.lookup<Scalar>("convergenceTolerance");
    if (config.convergenceTolerance <= S(0.0))
    {
        FatalError("SIMPLE.convergenceTolerance must be positive.");
    }

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
    config.turbulenceEnabled = turbulence.lookup<bool>("enabled");
    config.turbulenceModel = turbulence.lookup<std::string>("model");

    if
    (
        config.turbulenceEnabled
     && config.turbulenceModel != "kOmegaSST"
    )
    {
        FatalError
        (
            "Unsupported turbulence model: '"
          + config.turbulenceModel
          + "'. Only 'kOmegaSST' is supported."
        );
    }

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

    if (config.turbulenceEnabled)
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

    config.vtkOutputFilename = outputDict.lookup<std::string>("filename");
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

    readLinearSolvers(reader, config);

    return config;
}

} // namespace CaseConfig
