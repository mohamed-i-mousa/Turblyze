/******************************************************************************
 * @file CFDApplication.cpp
 * @brief Top-level application driver for the CFD solver
 *****************************************************************************/

#include "CFDApplication.h"

#include <iostream>
#include <array>
#include <map>
#include <set>
#include <algorithm>

#include <omp.h>
#include <eigen3/Eigen/Core>

#include "CaseReader.h"
#include "GradientScheme.h"
#include "SIMPLE.h"
#include "kOmegaSST.h"
#include "MeshReader.h"
#include "MeshChecker.h"
#include "VtkWriter.h"
#include "VtkBoundaryWriter.h"
#include "DerivedFields.h"
#include "LinearSolvers.h"
#include "Constraint.h"
#include "ErrorHandler.h"
#include "Logger.h"
#include <sstream>
#include <iomanip>

// ************************* Special Member Functions *************************

CFDApplication::CFDApplication(const std::string& caseFilePath)
:
    caseFilePath_(caseFilePath)
{}

CFDApplication::~CFDApplication() noexcept = default;

// *********************************** run ***********************************

void CFDApplication::run()
{
    loadCase();
    initParallelism();
    prepareMesh();
    setupBoundaryConditions();
    configureSolver();
    solve();
    postProcess();
    exportResults();
}

// ********************************* loadCase *********************************

void CFDApplication::loadCase()
{
    std::cout
        << '\n' << "--- 0. Loading Case ---" << '\n';

    caseReader_ = std::make_unique<CaseReader>(caseFilePath_);

    // Extract mesh configuration
    const auto& mesh = caseReader_->section("mesh");
    checkQuality_ = mesh.lookupOrDefault<bool>("checkQuality", true);

    // Extract parallelism settings
    if (caseReader_->hasSection("parallelism"))
    {
        const auto& parallelism = caseReader_->section("parallelism");
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
            numThreads_ = n;
        }
    }

    // Extract physical properties
    const auto& physicalProperties =
        caseReader_->section("physicalProperties");

    rho_ = physicalProperties.lookup<Scalar>("rho");
    mu_ = physicalProperties.lookup<Scalar>("mu");
    if (rho_ <= S(0)) FatalError("physicalProperties.rho must be positive.");
    if (mu_ <= S(0)) FatalError("physicalProperties.mu must be positive.");

    // Extract initial conditions
    const auto& initialConditions = caseReader_->section("initialConditions");
    initialVelocity_ = initialConditions.lookup<Vector>("U");
    initialPressure_ = initialConditions.lookup<Scalar>("p");

    debug_ =
        caseReader_->section("output").lookupOrDefault<bool>("debug", false);

    // Parse convection schemes
    convectionSchemes_ = parseConvectionSchemes();

    // Extract SIMPLE parameters
    const auto& simple = caseReader_->section("SIMPLE");

    maxIterations_ = simple.lookup<int>("numIterations");
    if (maxIterations_ <= 0)
    {
        FatalError("SIMPLE.numIterations must be a positive integer.");
    }

    convergenceTolerance_ = simple.lookup<Scalar>("convergenceTolerance");
    if (convergenceTolerance_ <= S(0))
    {
        FatalError("SIMPLE.convergenceTolerance must be positive.");
    }

    const auto& relaxFactors = simple.section("relaxationFactors");
    alphaU_ = relaxFactors.lookup<Scalar>("U");
    alphaP_ = relaxFactors.lookup<Scalar>("p");
    alphaK_ = relaxFactors.lookupOrDefault<Scalar>("k", S(0.5));
    alphaOmega_ = relaxFactors.lookupOrDefault<Scalar>("omega", S(0.5));

    if (alphaU_ <= S(0) || alphaU_ > S(2))
    {
        FatalError("SIMPLE.relaxationFactors.U must be in (0, 2].");
    }
    if (alphaP_ <= S(0) || alphaP_ > S(2))
    {
        FatalError("SIMPLE.relaxationFactors.p must be in (0, 2].");
    }
    if (alphaK_ <= S(0) || alphaK_ > S(2))
    {
        FatalError("SIMPLE.relaxationFactors.k must be in (0, 2].");
    }
    if (alphaOmega_ <= S(0) || alphaOmega_ > S(2))
    {
        FatalError("SIMPLE.relaxationFactors.omega must be in (0, 2].");
    }

    // Extract turbulence parameters
    const auto& turbulence = caseReader_->section("turbulence");
    turbulenceEnabled_ = turbulence.lookup<bool>("enabled");
    turbulenceModel_ = turbulence.lookup<std::string>("model");

    if (turbulenceEnabled_ && turbulenceModel_ != "kOmegaSST")
    {
        FatalError
        (
            "Unsupported turbulence model: '"
          + turbulenceModel_
          + "'. Only 'kOmegaSST' is supported."
        );
    }

    // Extract turbulence inlet/initial parameters
    turbIntensity_ =
        turbulence.lookupOrDefault<Scalar>("turbulenceIntensity", S(0.05));
    hydrDiameter_ =
        turbulence.lookupOrDefault<Scalar>("hydraulicDiameter", S(0.01));

    if (turbIntensity_ <= S(0))
    {
        FatalError("turbulence.turbulenceIntensity must be positive.");
    }
    if (hydrDiameter_ <= S(0))
    {
        FatalError("turbulence.hydraulicDiameter must be positive.");
    }

    if (turbulenceEnabled_)
    {
        // Module-computed defaults; user-supplied entries override them
        const Scalar defaultK =
            kOmegaSST::inletK(initialVelocity_, turbIntensity_);
        const Scalar defaultOmega =
            kOmegaSST::inletOmega(defaultK, hydrDiameter_);

        initialK_ =
            initialConditions.lookupOrDefault<Scalar>("k", defaultK);
        initialOmega_ =
            initialConditions.lookupOrDefault<Scalar>("omega", defaultOmega);

        if (initialK_ < S(0))
        {
            FatalError("initialConditions.k must be non-negative.");
        }
        if (initialOmega_ <= S(0))
        {
            FatalError("initialConditions.omega must be positive.");
        }
    }

    const auto& outputDict = caseReader_->section("output");
    vtkOutputFilename_ = outputDict.lookup<std::string>("filename");
    if (vtkOutputFilename_.empty())
    {
        FatalError("output.filename must not be empty.");
    }

    std::cout
        << "Case file loaded." << '\n';

    // Extract constraints (optional)
    if (caseReader_->hasSection("constraints"))
    {
        const auto& constraintsDict = caseReader_->section("constraints");

        if (constraintsDict.hasSection("velocity"))
        {
            const auto& velConstraint = constraintsDict.section("velocity");

            velocityConstraintEnabled_ = velConstraint.lookup<bool>("enabled");

            maxVelocityConstraint_ =
                velConstraint.lookup<Scalar>("maxVelocity");

            if (velocityConstraintEnabled_ && maxVelocityConstraint_ <= S(0))
            {
                FatalError
                (
                    "constraints.velocity.maxVelocity must be positive."
                );
            }
        }

        if (constraintsDict.hasSection("pressure"))
        {
            const auto& presConstraint = constraintsDict.section("pressure");

            pressureConstraintEnabled_ =
                presConstraint.lookup<bool>("enabled");

            minPressureConstraint_ =
                presConstraint.lookup<Scalar>("minPressure");

            maxPressureConstraint_ =
                presConstraint.lookup<Scalar>("maxPressure");

            if
            (
                pressureConstraintEnabled_
             && minPressureConstraint_ >= maxPressureConstraint_
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
}

// ****************************** initParallelism *****************************

void CFDApplication::initParallelism()
{
    omp_set_num_threads(numThreads_);
    Eigen::setNbThreads(numThreads_);

    std::cout
        << "OpenMP threads: " << numThreads_ << '\n';
}

// ******************************* prepareMesh *******************************

void CFDApplication::prepareMesh()
{
    std::cout
        << '\n' << "--- 1. Reading and Preparing Mesh ---" << '\n';

    const auto& mesh = caseReader_->section("mesh");
    const std::string meshFilePath = mesh.lookup<std::string>("file");

    MeshReader meshReader(meshFilePath);

    mesh_ = Mesh
    (
        meshReader.moveNodes(),
        meshReader.moveFaces(),
        meshReader.moveCells(),
        meshReader.moveBoundaryPatches()
    );

    std::cout
        << "Mesh Loaded: " << mesh_.numNodes() << " nodes, "
        << mesh_.numFaces() << " faces, " << mesh_.numCells()
        << " cells." << '\n';

    std::vector<FaceIntegrals> faceIntegrals(mesh_.numFaces());
    for (size_t faceIdx = 0; faceIdx < mesh_.numFaces(); ++faceIdx)
    {
        faceIntegrals[faceIdx] =
            mesh_.faces()[faceIdx].calculateGeometricProperties(mesh_.nodes());
    }
    if (debug_)
    {
        std::cout
            << "Geometric properties calculated for faces." << '\n';
    }

    {
        std::vector<Vector> approxCentroids(mesh_.numCells(), Vector{});
        std::vector<std::set<size_t>> cellNodes(mesh_.numCells());

        // Collect unique node indices for each cell
        for (const auto& face : mesh_.faces())
        {
            const size_t owner = face.ownerCell();
            for (size_t nodeIdx : face.nodeIndices())
            {
                cellNodes[owner].insert(nodeIdx);
            }

            if (!face.isBoundary())
            {
                const size_t neighbor = face.neighborCell().value();

                for (size_t nodeIdx : face.nodeIndices())
                {
                    cellNodes[neighbor].insert(nodeIdx);
                }
            }
        }

        // Compute centroid as average of node positions
        for (size_t cellIdx = 0; cellIdx < mesh_.numCells(); ++cellIdx)
        {
            if (!cellNodes[cellIdx].empty())
            {
                for (size_t nodeIdx : cellNodes[cellIdx])
                {
                    approxCentroids[cellIdx] += mesh_.nodes()[nodeIdx];
                }
                approxCentroids[cellIdx] /= S(cellNodes[cellIdx].size());
            }
        }

        // Check and correct any inverted face normals
        int flippedCount = 0;
        for (auto& face : mesh_.faces())
        {
            if (!face.isBoundary())
            {
                const Vector& ownerCell = approxCentroids[face.ownerCell()];

                const Vector& neighborCell =
                    approxCentroids[face.neighborCell().value()];

                const Vector dPN = neighborCell - ownerCell;

                if (dot(dPN, face.normal()) < 0)
                {
                    face.flipNormal();
                    flippedCount++;
                }
            }
        }

        if (debug_ && flippedCount > 0)
        {
            std::cout
                << "Corrected " << flippedCount
                << " inverted face normals." << '\n';
        }
    }

    for (auto& cell : mesh_.cells())
    {
        cell.calculateGeometricProperties(mesh_.faces(), faceIntegrals);
    }

    if (debug_)
    {
        std::cout
            << "Geometric properties calculated for cells." << '\n';
    }

    std::vector<Vector> cellCentroids(mesh_.numCells());

    for (size_t cellIdx = 0; cellIdx < mesh_.numCells(); ++cellIdx)
    {
        cellCentroids[cellIdx] = mesh_.cells()[cellIdx].centroid();
    }

    for (auto& face : mesh_.faces())
    {
        face.calculateDistanceProperties(cellCentroids);
    }

    if (debug_)
    {
        std::cout
            << "Distance properties calculated for faces." << '\n';
    }

    // Check mesh quality if requested
    if (checkQuality_)
    {
        MeshChecker meshChecker(mesh_);
        meshChecker.check();
    }
}

// ************************* setupBoundaryConditions *************************

void CFDApplication::setupBoundaryConditions()
{
    std::cout
        << '\n' << "--- 2. Setting Boundary Conditions ---" << '\n';

    for (const auto& patch : mesh_.patches())
    {
        bcManager_.addPatch(patch);
    }

    // Link boundary faces to their owning patches
    bcManager_.linkFaces(mesh_.faces());

    for (const auto& face : mesh_.faces())
    {
        if (face.isBoundary() && !face.patch().has_value())
        {
            FatalError
            (
                "Boundary face "
              + std::to_string(face.idx())
              + " has no patch after linking."
            );
        }
    }

    // Load boundary conditions from case file
    const auto& BCs = caseReader_->section("boundaryConditions");

    // Process velocity boundary conditions
    if (BCs.hasSection("U"))
    {
        const auto& velocityBCs = BCs.section("U");

        for (const auto& patchName : velocityBCs.sectionNames())
        {
            const auto& patchBC = velocityBCs.section(patchName);
            const std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                const Vector value = patchBC.lookup<Vector>("value");
                bcManager_.setFixedValue(patchName, Field::Ux, value.x());
                bcManager_.setFixedValue(patchName, Field::Uy, value.y());
                bcManager_.setFixedValue(patchName, Field::Uz, value.z());
            }
            else if (bcType == "noSlip")
            {
                bcManager_.setNoSlip(patchName, Field::Ux);
                bcManager_.setNoSlip(patchName, Field::Uy);
                bcManager_.setNoSlip(patchName, Field::Uz);
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, Field::Ux);
                bcManager_.setZeroGradient(patchName, Field::Uy);
                bcManager_.setZeroGradient(patchName, Field::Uz);
            }
            else
            {
                unknownBCType
                (
                    bcType, "U", patchName,
                    "fixedValue, noSlip, zeroGradient"
                );
            }
        }
    }

    // Process pressure boundary conditions
    bool hasFixedPressure = false;

    if (BCs.hasSection("p"))
    {
        const auto& pressureBCs = BCs.section("p");

        for (const auto& patchName : pressureBCs.sectionNames())
        {
            const auto& patchBC = pressureBCs.section(patchName);
            const std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                const Scalar value = patchBC.lookup<Scalar>("value");
                bcManager_.setFixedValue(patchName, Field::p, value);
                hasFixedPressure = true;
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, Field::p);
            }
            else
            {
                unknownBCType
                (
                    bcType, "p", patchName,
                    "fixedValue, zeroGradient"
                );
            }
        }

        // Derive pressure correction BCs from pressure BCs
        for (const auto& patchName : pressureBCs.sectionNames())
        {
            const auto& patchBC = pressureBCs.section(patchName);
            const std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                bcManager_.setFixedValue(patchName, Field::pCorr, S(0.0));
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, Field::pCorr);
            }
        }
    }

    if (!hasFixedPressure)
    {
        Warning
        (
            "No fixedValue pressure boundary condition found. "
            "The pressure field has no reference value, which "
            "may cause a singular pressure matrix."
        );
    }

    // Process turbulent kinetic energy BCs
    if (BCs.hasSection("k"))
    {
        const auto& kBCs = BCs.section("k");

        for (const auto& patchName : kBCs.sectionNames())
        {
            const auto& patchBC = kBCs.section(patchName);
            const std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                const std::string valStr =
                    patchBC.lookup<std::string>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    value = kOmegaSST::inletK(initialVelocity_, turbIntensity_);

                    std::cout
                        << "Inlet turbulence kinetic energy : " << value
                        << '\n';
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                bcManager_.setFixedValue
                (
                    patchName, Field::k, value
                );

            }
            else if (bcType == "kWallFunction")
            {
                bcManager_.setWallFunctionType
                (
                    patchName,
                    Field::k,
                    BCType::K_WALL_FUNCTION
                );
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, Field::k);
            }
            else
            {
                unknownBCType
                (
                    bcType, "k", patchName,
                    "fixedValue, kWallFunction, zeroGradient"
                );
            }
        }
    }

    // Process specific dissipation rate BCs
    if (BCs.hasSection("omega"))
    {
        const auto& omegaBCs = BCs.section("omega");

        for (const auto& patchName : omegaBCs.sectionNames())
        {
            const auto& patchBC = omegaBCs.section(patchName);
            const std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                const std::string valStr =
                    patchBC.lookup<std::string>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    // Use the patch k BC if fixed, else derive k from velocity
                    const BoundaryData& kPatchBC =
                        bcManager_.fieldBC(patchName, Field::k);

                    const Scalar kValue =
                        kPatchBC.type() == BCType::FIXED_VALUE
                      ? kPatchBC.fixedScalarValue()
                      : kOmegaSST::inletK(initialVelocity_, turbIntensity_);

                    value = kOmegaSST::inletOmega(kValue, hydrDiameter_);

                    std::cout
                        << "Inlet specific dissipation : " << value << '\n';
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                bcManager_.setFixedValue
                (
                    patchName, Field::omega, value
                );
            }
            else if (bcType == "omegaWallFunction")
            {
                bcManager_.setWallFunctionType
                (
                    patchName,
                    Field::omega,
                    BCType::OMEGA_WALL_FUNCTION
                );
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, Field::omega);
            }
            else
            {
                unknownBCType
                (
                    bcType, "omega", patchName,
                    "fixedValue, omegaWallFunction, zeroGradient"
                );
            }
        }
    }

    // Optional turbulent viscosity BC section
    if (BCs.hasSection("nut"))
    {
        const auto& nutBCs = BCs.section("nut");

        for (const auto& patchName : nutBCs.sectionNames())
        {
            const auto& patchBC = nutBCs.section(patchName);
            const std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                const Scalar value = patchBC.lookup<Scalar>("value");
                bcManager_.setFixedValue(patchName, Field::nut, value);
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, Field::nut);
            }
            else if (bcType == "nutWallFunction")
            {
                bcManager_.setWallFunctionType
                (
                    patchName,
                    Field::nut,
                    BCType::NUT_WALL_FUNCTION
                );
            }
            else
            {
                unknownBCType
                (
                    bcType, "nut", patchName,
                    "fixedValue, zeroGradient, nutWallFunction"
                );
            }
        }
    }

    // Validate patch names against mesh patches
    bcManager_.validatePatchNames();

    if (debug_)
    {
        bcManager_.printSummary();
    }

    std::cout
        << "Boundary conditions set for "
        << mesh_.patches().size() << " patches." << '\n';
}

// ***************************** configureSolver *****************************

void CFDApplication::configureSolver()
{
    gradScheme_ =
        std::make_unique<GradientScheme>(mesh_, bcManager_);

    solver_ =
        std::make_unique<SIMPLE>
        (
            mesh_,
            bcManager_,
            *gradScheme_,
            convectionSchemes_
        );

    solver_->setDebug(debug_);

    solver_->
        initialize
        (
            initialVelocity_,
            initialPressure_,
            initialK_,
            initialOmega_,
            turbulenceEnabled_
        );

    solver_->setPhysicalProperties(rho_, mu_);
    solver_->setRelaxationFactors(alphaU_, alphaP_, alphaK_, alphaOmega_);
    solver_->setConvergenceTolerance(convergenceTolerance_);
    solver_->setMaxIterations(maxIterations_);

    // Linear solver settings (defaults per equation)
    std::string uSolverName     = "BiCGSTAB";
    std::string pSolverName     = "PCG";
    std::string kSolverName     = "BiCGSTAB";
    std::string omegaSolverName = "BiCGSTAB";
    std::string uPreconditioner     = "Jacobi";
    std::string pPreconditioner     = "Jacobi";
    std::string kPreconditioner     = "Jacobi";
    std::string omegaPreconditioner = "Jacobi";
    Scalar uTol     = S(1e-8);
    Scalar pTol     = S(1e-6);
    Scalar kTol     = S(1e-6);
    Scalar omegaTol = S(1e-6);
    int uMaxIter     = 1000;
    int pMaxIter     = 1000;
    int kMaxIter     = 1000;
    int omegaMaxIter = 1000;

    if (caseReader_->hasSection("linearSolvers"))
    {
        const auto& solvers = caseReader_->section("linearSolvers");

        readAndValidateSolverConfig
        (
            solvers,
            "U",
            uSolverName,
            uPreconditioner,
            uTol,
            uMaxIter
        );
        readAndValidateSolverConfig
        (
            solvers,
            "p",
            pSolverName,
            pPreconditioner,
            pTol,
            pMaxIter
        );

        if (turbulenceEnabled_)
        {
            readAndValidateSolverConfig
            (
                solvers,
                "k",
                kSolverName,
                kPreconditioner,
                kTol,
                kMaxIter
            );
            readAndValidateSolverConfig
            (
                solvers,
                "omega",
                omegaSolverName,
                omegaPreconditioner,
                omegaTol,
                omegaMaxIter
            );
        }
    }

    solver_->setMomentumSolver
    (
        createLinearSolver(uSolverName, uTol, uMaxIter)
    );
    solver_->setPressureSolver
    (
        createLinearSolver(pSolverName, pTol, pMaxIter)
    );

    if (turbulenceEnabled_)
    {
        solver_->setTurbulenceSolvers
        (
            createLinearSolver(kSolverName, kTol, kMaxIter),
            createLinearSolver(omegaSolverName, omegaTol, omegaMaxIter)
        );
    }

    Constraint& constraintSystem = solver_->constraintSystem();

    if (velocityConstraintEnabled_)
    {
        constraintSystem.setVelocityConstraints(maxVelocityConstraint_);
    }
    if (pressureConstraintEnabled_)
    {
        constraintSystem.setPressureConstraints
        (
            minPressureConstraint_,
            maxPressureConstraint_
        );
    }
    constraintSystem.enableConstraints
    (
        velocityConstraintEnabled_,
        pressureConstraintEnabled_
    );

    // Framed Phase 3 summary — always printed, regardless of debug mode
    Logger::sectionHeader("Initializing SIMPLE Solver");

    Logger::subsection("Linear solvers");
    Logger::linearSolverConfigHeader();
    Logger::linearSolverConfigRow("U", uSolverName, uTol, uMaxIter);
    Logger::linearSolverConfigRow("p", pSolverName, pTol, pMaxIter);
    if (turbulenceEnabled_)
    {
        Logger::linearSolverConfigRow("k", kSolverName, kTol, kMaxIter);
        Logger::linearSolverConfigRow
        (
            "omega",
            omegaSolverName,
            omegaTol,
            omegaMaxIter
        );
    }

    Logger::subsection("SIMPLE controls");
    Logger::keyValue("Max iterations", maxIterations_);
    Logger::keyValue("Convergence tolerance", convergenceTolerance_);
    Logger::keyValue("Velocity relaxation", alphaU_);
    Logger::keyValue("Pressure relaxation", alphaP_);
    if (turbulenceEnabled_)
    {
        Logger::keyValue("k relaxation", alphaK_);
        Logger::keyValue("omega relaxation", alphaOmega_);
    }

    if (velocityConstraintEnabled_ || pressureConstraintEnabled_)
    {
        Logger::subsection("Field constraints");
        if (velocityConstraintEnabled_)
        {
            std::ostringstream os;
            os << std::scientific << std::setprecision(6)
               << maxVelocityConstraint_ << " m/s";
            Logger::keyValue("Velocity max", os.str());
        }
        if (pressureConstraintEnabled_)
        {
            std::ostringstream os;
            os << "[" << std::scientific << std::setprecision(6)
               << minPressureConstraint_ << ", "
               << maxPressureConstraint_ << "] Pa";
            Logger::keyValue("Pressure range", os.str());
        }
    }

    if (turbulenceEnabled_)
    {
        Logger::subsection("Turbulence initialization");
        Logger::keyValue("Model", turbulenceModel_);
        Logger::keyValue
        (
            "Wall distance",
            std::string
            {
                solver_->wallDistanceConverged()
              ? "meshWave converged"
              : "meshWave hit iteration cap (results may be degraded)"
            }
        );
        Logger::keyValue("Fields initialized", std::string{"k, omega, nut"});
    }

    Logger::iterationFooter();
}

// ********************************** solve ***********************************

void CFDApplication::solve()
{
    std::cout
        << '\n' << "--- 4. Solving Steady-State Flow with SIMPLE ---"
        << '\n';

    solver_->solve();
}

// ******************************** postProcess *******************************

void CFDApplication::postProcess()
{
    std::cout
         << '\n' << "--- 5. Extracting Solution Fields ---" << '\n';

    const ScalarField& Ux = solver_->Ux();
    const ScalarField& Uy = solver_->Uy();
    const ScalarField& Uz = solver_->Uz();
    const ScalarField& pressure = solver_->pressure();

    if (debug_)
    {
        std::cout
            << "Solution extracted." << '\n';
    }

    std::cout
         << '\n' << "--- 6. Post-Processing Results ---" << '\n';

    // Calculate velocity magnitude
    const ScalarField velocityMag = VTK::velocityMagnitude(Ux, Uy, Uz);

    if (Ux.size() == 0)
    {
        Warning("Solution fields are empty. Skipping statistics.");
        return;
    }

    // Print statistics
    Scalar maximumVelocity = S(0.0);
    Scalar averageVelocity = S(0.0);
    Scalar maximumPressure = pressure[0];
    Scalar minimumPressure = pressure[0];

    for (size_t cellIdx = 0; cellIdx < Ux.size(); ++cellIdx)
    {
        const Scalar vmag = velocityMag[cellIdx];
        maximumVelocity = std::max(maximumVelocity, vmag);
        averageVelocity += vmag;

        maximumPressure = std::max(maximumPressure, pressure[cellIdx]);
        minimumPressure = std::min(minimumPressure, pressure[cellIdx]);
    }
    averageVelocity /= S(Ux.size());

    std::cout
        << "Flow Statistics:" << '\n';
    std::cout
        << "  Max velocity magnitude: " << maximumVelocity
        << " m/s" << '\n';
    std::cout
        << "  Average velocity magnitude: "
        << averageVelocity << " m/s" << '\n';
    std::cout
        << "  Pressure range: [" << minimumPressure
        << ", " << maximumPressure << "] Pa" << '\n';
}

// ****************************** exportResults *******************************

void CFDApplication::exportResults()
{
    std::cout
         << '\n' << "--- 7. Exporting Results to VTK ---" << '\n';

    const ScalarField& Ux = solver_->Ux();
    const ScalarField& Uy = solver_->Uy();
    const ScalarField& Uz = solver_->Uz();
    const ScalarField& pressure = solver_->pressure();

    // Calculate velocity magnitude
    const ScalarField velocityMag = VTK::velocityMagnitude(Ux, Uy, Uz);

    // Prepare scalar fields for export
    std::map<std::string, const ScalarField*>
    scalarFieldsToVtk;

    scalarFieldsToVtk["pressure"] = &pressure;

    scalarFieldsToVtk["velocityMagnitude"] = &velocityMag;

    // Add turbulence fields if available
    if (turbulenceEnabled_)
    {
        scalarFieldsToVtk["k"] = &solver_->turbulentKineticEnergy();
        scalarFieldsToVtk["omega"] = &solver_->specificDissipationRate();
        scalarFieldsToVtk["nut"] = &solver_->turbulentViscosity();
        scalarFieldsToVtk["wallDistance"] = &solver_->wallDistance();
    }

    // Prepare vector fields for export
    std::map<std::string, std::array<const ScalarField*, 3>>
    vectorFieldsToVtk;

    vectorFieldsToVtk["velocity"] = {&Ux, &Uy, &Uz};

    // Determine VTU filename
    std::string vtuFilename = vtkOutputFilename_;
    const std::string ext = ".vtu";
    if (!vtuFilename.ends_with(ext))
    {
        vtuFilename += ext;
    }

    if (debug_)
    {
        std::cout
            << '\n' << "Exporting results to VTK UnstructuredGrid..." << '\n';
    }

    VTK::writeVtkUnstructuredGrid
    (
        vtuFilename,
        mesh_,
        scalarFieldsToVtk,
        vectorFieldsToVtk,
        debug_
    );

    // Export wall boundary data (yPlus, wallShearStress) as VTP
    if (turbulenceEnabled_)
    {
        std::string vtpFilename = vtuFilename;
        const size_t dotPos = vtpFilename.rfind(".vtu");

        if (dotPos != std::string::npos)
        {
            vtpFilename.replace(dotPos, 4, "_wall.vtp");
        }
        else
        {
            vtpFilename += "_wall.vtp";
        }

        std::map<std::string, const FaceData<Scalar>*>
        wallScalarFields;

        wallScalarFields["yPlus"] = &solver_->yPlus();
        wallScalarFields["wallShearStress"] = &solver_->wallShearStress();

        VTK::writeWallBoundaryData
        (
            vtpFilename,
            mesh_,
            wallScalarFields,
            debug_
        );
    }

    std::cout
        << '\n' << "=== CFD Results Exported Successfully ===" << '\n';
    std::cout
        << "File: " << vtuFilename << '\n';
}

// ************************* createConvectionScheme ***********************

std::unique_ptr<ConvectionSchemes>
CFDApplication::createConvectionScheme(const std::string& name)
{
    if (name == "Upwind")
    {
        return std::make_unique<UpwindScheme>();
    }

    if (name == "CentralDifference")
    {
        return std::make_unique<CentralDifferenceScheme>();
    }

    if (name == "SecondOrderUpwind")
    {
        return std::make_unique<SecondOrderUpwindScheme>();
    }

    FatalError("Unknown convection scheme: " + name);
}

// ************************** createLinearSolver **************************

std::unique_ptr<LinearSolver> CFDApplication::createLinearSolver
(
    std::string_view name,
    Scalar tolerance,
    int maxIterations
)
{
    if (name == BiCGSTAB::typeName)
    {
        return std::make_unique<BiCGSTAB>(tolerance, maxIterations);
    }

    if (name == PCG::typeName)
    {
        return std::make_unique<PCG>(tolerance, maxIterations);
    }

    FatalError
    (
        "Unknown linear solver '" + std::string(name)
      + "'. Supported solvers: BiCGSTAB, PCG."
    );
}

// ***************************** unknownBCType *****************************

void CFDApplication::unknownBCType
(
    std::string_view bcType,
    std::string_view fieldName,
    std::string_view patchName,
    std::string_view validList
)
{
    FatalError
    (
        "Unknown boundary condition type '" + std::string(bcType)
      + "' for field '" + std::string(fieldName)
      + "' on patch '" + std::string(patchName)
      + "'. Valid types: " + std::string(validList)
    );
}

// ********************* readAndValidateSolverConfig **********************

void CFDApplication::readAndValidateSolverConfig
(
    const CaseReader& solvers,
    const std::string& key,
    std::string& solverName,
    std::string& preconditioner,
    Scalar& tolerance,
    int& maxIterations
)
{
    if (solvers.hasSection(key))
    {
        const auto& s = solvers.section(key);
        solverName = s.lookupOrDefault<std::string>("solver", solverName);
        preconditioner =
            s.lookupOrDefault<std::string>
            (
                "preconditioner",
                preconditioner
            );
        tolerance =
            s.lookupOrDefault<Scalar>("tolerance", tolerance);
        maxIterations =
            s.lookupOrDefault<int>("maxIter", maxIterations);
    }

    if (tolerance <= S(0))
    {
        FatalError("linearSolvers." + key + ".tolerance must be positive.");
    }
    if (maxIterations <= 0)
    {
        FatalError
        (
            "linearSolvers." + key + ".maxIter must be a positive integer."
        );
    }
}

// ************************ parseConvectionSchemes ************************

ConvectionScheme CFDApplication::parseConvectionSchemes() const
{
    const auto& schemesDict = caseReader_->section("numericalSchemes");

    ConvectionScheme schemes;

    if (schemesDict.hasSection("convection"))
    {
        const auto& convSec = schemesDict.section("convection");

        const std::string defaultSchemeName =
            convSec.lookupOrDefault<std::string>("default", "Upwind");

        schemes.defaultScheme = createConvectionScheme(defaultSchemeName);

        if (debug_)
        {
            std::cout
                << "Default convection scheme: "
                << defaultSchemeName << '\n';
        }

        // Per-equation overrides (fall back to default if unset)
        const std::string uSchemeName =
            convSec.lookupOrDefault<std::string>("U", "");

        if (!uSchemeName.empty())
        {
            schemes.momentumScheme =
                createConvectionScheme(uSchemeName);

            if (debug_)
            {
                std::cout
                    << "Momentum convection scheme: "
                    << uSchemeName << '\n';
            }
        }

        const std::string kSchemeName =
            convSec.lookupOrDefault<std::string>("k", "");

        if (!kSchemeName.empty())
        {
            schemes.kScheme =
                createConvectionScheme(kSchemeName);

            if (debug_)
            {
                std::cout
                    << "k convection scheme: "
                    << kSchemeName << '\n';
            }
        }

        const std::string omegaSchemeName =
            convSec.lookupOrDefault<std::string>("omega", "");

        if (!omegaSchemeName.empty())
        {
            schemes.omegaScheme =
                createConvectionScheme(omegaSchemeName);

            if (debug_)
            {
                std::cout
                    << "omega convection scheme: "
                    << omegaSchemeName << '\n';
            }
        }
    }
    else
    {
        FatalError
        (
            "Missing 'convection' sub-section "
            "in numericalSchemes"
        );
    }

    return schemes;
}
