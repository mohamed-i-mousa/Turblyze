/******************************************************************************
 * @file CFDApplication.cpp
 * @brief Top-level application driver for the CFD solver
 *****************************************************************************/

#include "CFDApplication.h"

#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

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


// ******************************* Constructor *******************************

CFDApplication::CFDApplication(const std::string& caseFilePath)
:
    caseFilePath_(caseFilePath)
{}

CFDApplication::~CFDApplication() = default;


// *********************************** run ***********************************

void CFDApplication::run()
{
    loadCase();
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
        << std::endl << "--- 0. Loading Case ---" << std::endl;

    caseReader_ = std::make_unique<CaseReader>(caseFilePath_);

    // Extract mesh configuration
    const auto& mesh = caseReader_->section("mesh");
    checkQuality_ = mesh.lookupOrDefault<bool>("checkQuality", true);

    // Extract physical properties
    const auto& physicalProperties =
        caseReader_->section("physicalProperties");

    rho_ = physicalProperties.lookup<Scalar>("rho");
    mu_ = physicalProperties.lookup<Scalar>("mu");

    // Extract initial conditions
    const auto& initialConditions = caseReader_->section("initialConditions");
    initialVelocity_ = initialConditions.lookup<Vector>("U");
    initialPressure_ = initialConditions.lookup<Scalar>("p");

    // Parse convection schemes
    convectionSchemes_ = parseConvectionSchemes();

    // Extract SIMPLE parameters
    const auto& simple = caseReader_->section("SIMPLE");

    maxIterations_ = simple.lookup<int>("numIterations");

    convergenceTolerance_ = simple.lookup<Scalar>("convergenceTolerance");

    const auto& relaxFactors = simple.section("relaxationFactors");
    alphaU_ = relaxFactors.lookup<Scalar>("U");
    alphaP_ = relaxFactors.lookup<Scalar>("p");
    alphaK_ = relaxFactors.lookupOrDefault<Scalar>("k", 0.5);
    alphaOmega_ = relaxFactors.lookupOrDefault<Scalar>("omega", 0.5);

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

    // Compute turbulence initial conditions
    turbIntensity_ =
        turbulence.lookupOrDefault<Scalar>("turbulenceIntensity", S(0.05));
    hydrDiameter_ =
        turbulence.lookupOrDefault<Scalar>("hydraulicDiameter", S(0.01));

    Scalar lTurb = std::max(S(0.07) * hydrDiameter_, smallValue);
    Scalar UMag = initialVelocity_.magnitude();

    defaultK_ =
        std::max
        (
            S(1.5) * (turbIntensity_ * UMag) * (turbIntensity_ * UMag),
            S(1e-8)
        );

    defaultOmega_ =
        std::max
        (
            std::sqrt(defaultK_)
          / (std::pow(kOmegaSST::coeffs_.betaStar, S(0.25)) * lTurb),
            S(1e-4)
        );

    initialK_ = initialConditions.lookupOrDefault<Scalar>("k", defaultK_);

    initialOmega_ =
        initialConditions.lookupOrDefault<Scalar>("omega", defaultOmega_);

    // Extract output configuration
    const auto& outputDict = caseReader_->section("output");
    vtkOutputFilename_ = outputDict.lookup<std::string>("filename");
    debug_ = outputDict.lookupOrDefault<bool>("debug", false);

    std::cout
        << "Case file loaded." << std::endl;

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
        }
    }
}


// ******************************* prepareMesh *******************************

void CFDApplication::prepareMesh()
{
    std::cout
        << std::endl << "--- 1. Reading and Preparing Mesh ---" << std::endl;

    const auto& mesh = caseReader_->section("mesh");
    std::string meshFilePath = mesh.lookup<std::string>("file");

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
        << " cells." << std::endl;

    std::vector<FaceIntegrals> faceIntegrals(mesh_.numFaces());
    for (size_t faceIdx = 0; faceIdx < mesh_.numFaces(); ++faceIdx)
    {
        faceIntegrals[faceIdx] =
            mesh_.faces()[faceIdx].calculateGeometricProperties(mesh_.nodes());
    }
    if (debug_)
    {
        std::cout
            << "Geometric properties calculated for faces."
            << std::endl;
    }

    {
        std::vector<Vector> approxCentroids(mesh_.numCells(), Vector{});
        std::vector<std::set<size_t>> cellNodes(mesh_.numCells());

        // Collect unique node indices for each cell
        for (const auto& face : mesh_.faces())
        {
            size_t owner = face.ownerCell();
            for (size_t nodeIdx : face.nodeIndices())
            {
                cellNodes[owner].insert(nodeIdx);
            }

            if (!face.isBoundary())
            {
                size_t neighbor = face.neighborCell().value();

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

                Vector dPN = neighborCell - ownerCell;

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
                << " inverted face normals." << std::endl;
        }
    }

    for (auto& cell : mesh_.cells())
    {
        cell.calculateGeometricProperties(mesh_.faces(), faceIntegrals);
    }

    if (debug_)
    {
        std::cout
            << "Geometric properties calculated for cells."
            << std::endl;
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
            << "Distance properties calculated for faces."
            << std::endl;
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
        << std::endl << "--- 2. Setting Boundary Conditions ---" << std::endl;

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
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                Vector value = patchBC.lookup<Vector>("value");
                bcManager_.setFixedValue(patchName, "U", value);
            }
            else if (bcType == "noSlip")
            {
                bcManager_.setNoSlip(patchName, "U");
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "U");
            }
            else
            {
                FatalError
                (
                    "Unknown boundary condition type '"
                  + bcType
                  + "' for field 'U' on patch '"
                  + patchName
                  + "'. Valid types: "
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
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                Scalar value = patchBC.lookup<Scalar>("value");
                bcManager_.setFixedValue(patchName, "p", value);
                hasFixedPressure = true;
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "p");
            }
            else
            {
                FatalError
                (
                    "Unknown boundary condition type '"
                  + bcType
                  + "' for field 'p' on patch '"
                  + patchName
                  + "'. Valid types: "
                    "fixedValue, zeroGradient"
                );
            }
        }

        // Derive pressure correction BCs from pressure BCs
        for (const auto& patchName : pressureBCs.sectionNames())
        {
            const auto& patchBC = pressureBCs.section(patchName);
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                bcManager_.setFixedValue(patchName, "pCorr", S(0.0));
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "pCorr");
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
    const Scalar lengthScale = std::max(S(0.07) * hydrDiameter_, S(1e-20));

    if (BCs.hasSection("k"))
    {
        const auto& kBCs = BCs.section("k");

        for (const auto& patchName : kBCs.sectionNames())
        {
            const auto& patchBC = kBCs.section(patchName);
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                std::string valStr = patchBC.lookup<std::string>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    const BoundaryData& velocityBC =
                        bcManager_.fieldBC(patchName, "U");

                    Scalar UMag = initialVelocity_.magnitude();

                    if
                    (
                        velocityBC.type() == BCType::FIXED_VALUE
                     && velocityBC.valueType() == BCValueType::VECTOR
                    )
                    {
                        UMag = velocityBC.fixedVectorValue().magnitude();
                    }
                    else if (velocityBC.type() == BCType::NO_SLIP)
                    {
                        UMag = S(0.0);
                    }

                    const Scalar uPrime = turbIntensity_ * UMag;
                    value = std::max(S(1.5) * uPrime * uPrime, S(1e-8));

                    std::cout
                        << "Inlet turbulence kinetic energy : " << value
                        << std::endl;
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                bcManager_.setFixedValue
                (
                    patchName, "k", value
                );

            }
            else if (bcType == "kWallFunction")
            {
                bcManager_.setWallFunctionType
                (
                    patchName,
                    "k",
                    BCType::K_WALL_FUNCTION
                );
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "k");
            }
            else
            {
                FatalError
                (
                    "Unknown boundary condition type '"
                  + bcType
                  + "' for field 'k' on patch '"
                  + patchName
                  + "'. Valid types: "
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
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                std::string valStr = patchBC.lookup<std::string>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    // Determine k value for this patch
                    Scalar kValue = S(0.0);

                    const BoundaryData& kPatchBC =
                        bcManager_.fieldBC(patchName, "k");

                    if
                    (
                        kPatchBC.type() == BCType::FIXED_VALUE
                     && kPatchBC.valueType() == BCValueType::SCALAR
                    )
                    {
                        kValue = kPatchBC.fixedScalarValue();
                    }
                    else
                    {
                        // Compute k from velocity BC
                        const BoundaryData& velocityBC =
                            bcManager_.fieldBC(patchName, "U");

                        Scalar UMag = initialVelocity_.magnitude();

                        if
                        (
                            velocityBC.type() == BCType::FIXED_VALUE
                         && velocityBC.valueType() == BCValueType::VECTOR
                        )
                        {
                            UMag = velocityBC.fixedVectorValue().magnitude();
                        }
                        else if (velocityBC.type() == BCType::NO_SLIP)
                        {
                            UMag = S(0.0);
                        }

                        const Scalar uPrime = turbIntensity_ * UMag;
                        kValue = std::max(S(1.5) * uPrime * uPrime, S(1e-8));
                    }

                    // Compute omega from k
                    const Scalar omegaValue =
                        std::sqrt(std::max(kValue, S(0.0)))
                      / (std::pow(kOmegaSST::coeffs_.betaStar, S(0.25))
                       * lengthScale);

                    value = std::max(omegaValue, S(1e-4));

                    std::cout
                        << "Inlet specific dissipation : " << value
                        << std::endl;
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                bcManager_.setFixedValue
                (
                    patchName, "omega", value
                );
            }
            else if (bcType == "omegaWallFunction")
            {
                bcManager_.setWallFunctionType
                (
                    patchName,
                    "omega",
                    BCType::OMEGA_WALL_FUNCTION
                );
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "omega");
            }
            else
            {
                FatalError
                (
                    "Unknown boundary condition type '"
                  + bcType
                  + "' for field 'omega' on patch '"
                  + patchName
                  + "'. Valid types: "
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
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                Scalar value = patchBC.lookup<Scalar>("value");
                bcManager_.setFixedValue(patchName, "nut", value);
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "nut");
            }
            else if (bcType == "nutWallFunction")
            {
                bcManager_.setWallFunctionType
                (
                    patchName,
                    "nut",
                    BCType::NUT_WALL_FUNCTION
                );
            }
            else
            {
                FatalError
                (
                    "Unknown boundary condition type '"
                  + bcType
                  + "' for field 'nut' on patch '"
                  + patchName
                  + "'. Valid types: "
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
        << mesh_.patches().size() << " patches." << std::endl;
}


// ***************************** configureSolver *****************************

void CFDApplication::configureSolver()
{
    std::cout
        << std::endl << "--- 3. Initializing SIMPLE Solver ---" << std::endl;

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

    // Set debug mode before any output-producing calls
    solver_->setDebug(debug_);

    if (debug_ && turbulenceEnabled_)
    {
        std::cout
            << "Turbulence modeling enabled: "
            << turbulenceModel_ << std::endl;
    }

    // Initialize all solution fields
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

    // Configure SIMPLE parameters
    solver_->setRelaxationFactors(alphaU_, alphaP_, alphaK_, alphaOmega_);
    solver_->setConvergenceTolerance(convergenceTolerance_);
    solver_->setMaxIterations(maxIterations_);

    // Configure linear solvers from case file
    if (caseReader_->hasSection("linearSolvers"))
    {
        const auto& solvers = caseReader_->section("linearSolvers");

        // Momentum solver (U)
        if (solvers.hasSection("U"))
        {
            const auto& USection = solvers.section("U");
            LinearSolver momentumSolver
            (
                "momentum",
                USection.lookupOrDefault<Scalar>("tolerance", S(1e-8)),
                USection.lookupOrDefault<int>("maxIter", 1000)
            );
            momentumSolver.setILUTParameters
            (
                USection.lookupOrDefault<int>("fillFactor", 5),
                USection.lookupOrDefault<Scalar>("dropTol", S(1e-3))
            );
            solver_->setMomentumSolver(momentumSolver);

            if (debug_)
            {
                std::cout
                    << "  Momentum solver: tolerance="
                    << USection.lookupOrDefault<Scalar>(
                        "tolerance", S(1e-8))
                    << ", maxIter="
                    << USection.lookupOrDefault<int>(
                        "maxIter", 1000)
                    << std::endl;
            }
        }

        // Pressure solver (p)
        if (solvers.hasSection("p"))
        {
            const auto& pSection = solvers.section("p");
            LinearSolver pressureSolver
            (
                "pCorr",
                pSection.lookupOrDefault<Scalar>("tolerance", S(1e-6)),
                pSection.lookupOrDefault<int>("maxIter", 1000)
            );

            Scalar initialShift =
                pSection.lookupOrDefault<Scalar>("initialShift", S(1e-2));

            pressureSolver.setICParameters(initialShift);
            solver_->setPressureSolver(pressureSolver);

            if (debug_)
            {
                std::cout
                    << "  Pressure solver: tolerance="
                    << pSection.lookupOrDefault<Scalar>(
                        "tolerance", S(1e-6))
                    << ", maxIter="
                    << pSection.lookupOrDefault<int>(
                        "maxIter", 1000)
                    << std::endl;
            }
        }

        // Turbulence solvers (k, omega)
        if (turbulenceEnabled_)
        {
            LinearSolver kSolver
            (
                "k",
                solvers.hasSection("k")
              ? solvers.section("k").lookupOrDefault<Scalar>
                ("tolerance", S(1e-6))
              : S(1e-6),

                solvers.hasSection("k")
              ? solvers.section("k").lookupOrDefault<int>
                ("maxIter", 1000)
              : 1000
            );
            kSolver.setILUTParameters
            (
                solvers.hasSection("k")
              ? solvers.section("k").lookupOrDefault<int>
                ("fillFactor", 5)
              : 5,

                solvers.hasSection("k")
              ? solvers.section("k").lookupOrDefault<Scalar>
                ("dropTol", S(1e-3))
              : S(1e-3)
            );

            LinearSolver omegaSolver
            (
                "omega",
                solvers.hasSection("omega")
              ? solvers.section("omega").lookupOrDefault<Scalar>
                ("tolerance", S(1e-6))
              : S(1e-6),

                solvers.hasSection("omega")
              ? solvers.section("omega").lookupOrDefault<int>
                ("maxIter", 1000)
              : 1000
            );
            omegaSolver.setILUTParameters
            (
                solvers.hasSection("omega")
              ? solvers.section("omega").lookupOrDefault<int>
                ("fillFactor", 5)
              : 5,

                solvers.hasSection("omega")
              ? solvers.section("omega").lookupOrDefault<Scalar>
                ("dropTol", S(1e-3))
              : S(1e-3)
            );

            solver_->setTurbulenceSolvers(kSolver, omegaSolver);
        }
    }

    if (debug_)
    {
        std::cout
            << "SIMPLE parameters:" << std::endl;
        std::cout
            << "  Max iterations: " << maxIterations_
            << std::endl;
        std::cout
            << "  Convergence tolerance: "
            << convergenceTolerance_ << std::endl;
        std::cout
            << "  Velocity relaxation: " << alphaU_
            << std::endl;
        std::cout
            << "  Pressure relaxation: " << alphaP_
            << std::endl;

        if (turbulenceEnabled_)
        {
            std::cout
                << "  k relaxation: " << alphaK_
                << std::endl;
            std::cout
                << "  omega relaxation: " << alphaOmega_
                << std::endl;
        }
    }

    // Configure field constraints
    Constraint& constraintSystem = solver_->constraintSystem();

    if (velocityConstraintEnabled_)
    {
        constraintSystem.setVelocityConstraints(
            maxVelocityConstraint_);

        if (debug_)
        {
            std::cout
                << "Velocity constraint enabled: max = "
                << maxVelocityConstraint_ << " m/s"
                << std::endl;
        }
    }
    if (pressureConstraintEnabled_)
    {
        constraintSystem.setPressureConstraints
        (
            minPressureConstraint_,
            maxPressureConstraint_
        );

        if (debug_)
        {
            std::cout
                << "Pressure constraint enabled: ["
                << minPressureConstraint_
                << ", " << maxPressureConstraint_
                << "] Pa" << std::endl;
        }
    }
    constraintSystem.enableConstraints
    (
        velocityConstraintEnabled_,
        pressureConstraintEnabled_
    );

    std::cout
        << "SIMPLE solver initialized." << std::endl;
}


// ********************************** solve ***********************************

void CFDApplication::solve()
{
    std::cout
        << std::endl << "--- 4. Solving Steady-State Flow with SIMPLE ---"
        << std::endl;

    solver_->solve();
}


// ******************************** postProcess *******************************

void CFDApplication::postProcess()
{
    std::cout
         << std::endl << "--- 5. Extracting Solution Fields ---" << std::endl;

    const VectorField& velocity = solver_->velocity();
    const ScalarField& pressure = solver_->pressure();

    if (debug_)
    {
        std::cout
            << "Solution extracted." << std::endl;
    }

    std::cout
         << std::endl << "--- 6. Post-Processing Results ---" << std::endl;

    // Calculate velocity magnitude
    ScalarField velocityMag = VTK::velocityMagnitude(velocity);

    if (velocity.size() == 0)
    {
        Warning("Solution fields are empty. Skipping statistics.");
        return;
    }

    // Print statistics
    Scalar maximumVelocity = 0.0;
    Scalar averageVelocity = 0.0;
    Scalar maximumPressure = pressure[0];
    Scalar minimumPressure = pressure[0];

    for (size_t cellIdx = 0; cellIdx < velocity.size(); ++cellIdx)
    {
        Scalar vmag = velocityMag[cellIdx];
        maximumVelocity = std::max(maximumVelocity, vmag);
        averageVelocity += vmag;

        maximumPressure = std::max(maximumPressure, pressure[cellIdx]);
        minimumPressure = std::min(minimumPressure, pressure[cellIdx]);
    }
    averageVelocity /= S(velocity.size());

    std::cout
        << "Flow Statistics:" << std::endl;
    std::cout
        << "  Max velocity magnitude: " << maximumVelocity
        << " m/s" << std::endl;
    std::cout
        << "  Average velocity magnitude: "
        << averageVelocity << " m/s" << std::endl;
    std::cout
        << "  Pressure range: [" << minimumPressure
        << ", " << maximumPressure << "] Pa" << std::endl;
}


// ****************************** exportResults *******************************

void CFDApplication::exportResults()
{
    std::cout
         << std::endl << "--- 7. Exporting Results to VTK ---" << std::endl;

    const VectorField& velocity = solver_->velocity();
    const ScalarField& pressure = solver_->pressure();

    // Calculate velocity magnitude
    ScalarField velocityMag = VTK::velocityMagnitude(velocity);

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
    std::map<std::string, const VectorField*>
    vectorFieldsToVtk;

    vectorFieldsToVtk["velocity"] = &velocity;

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
            << std::endl
            << "Exporting results to VTK UnstructuredGrid..."
            << std::endl;
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
        size_t dotPos = vtpFilename.rfind(".vtu");

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

        wallScalarFields["yPlus"] =
            &solver_->yPlus();
        wallScalarFields["wallShearStress"] =
            &solver_->wallShearStress();

        VTK::writeWallBoundaryData
        (
            vtpFilename,
            mesh_,
            wallScalarFields,
            debug_
        );
    }

    std::cout
        << std::endl << "=== CFD Results Exported Successfully ==="
        << std::endl;
    std::cout
        << "File: " << vtuFilename << std::endl;
}


// ************************* createConvectionScheme ***********************

std::unique_ptr<ConvectionSchemes>
CFDApplication::createConvectionScheme(const std::string& name)
{
    if (name == "Upwind")
    {
        return std::make_unique<UpwindScheme>();
    }
    else if (name == "CentralDifference")
    {
        return std::make_unique<CentralDifferenceScheme>();
    }
    else if (name == "SecondOrderUpwind")
    {
        return std::make_unique<SecondOrderUpwindScheme>();
    }
    else
    {
        FatalError("Unknown convection scheme: " + name);
    }
}


// ************************ parseConvectionSchemes ************************

ConvectionScheme CFDApplication::parseConvectionSchemes()
{
    const auto& schemesDict = caseReader_->section("numericalSchemes");

    ConvectionScheme schemes;

    if (schemesDict.hasSection("convection"))
    {
        const auto& convSec = schemesDict.section("convection");

        std::string defaultSchemeName =
            convSec.lookupOrDefault<std::string>("default", "Upwind");

        schemes.defaultScheme = createConvectionScheme(defaultSchemeName);

        if (debug_)
        {
            std::cout
                << "Default convection scheme: "
                << defaultSchemeName << std::endl;
        }

        // Per-equation overrides (fall back to default if unset)
        std::string uSchemeName =
            convSec.lookupOrDefault<std::string>("U", "");

        if (!uSchemeName.empty())
        {
            schemes.momentumScheme =
                createConvectionScheme(uSchemeName);

            if (debug_)
            {
                std::cout
                    << "Momentum convection scheme: "
                    << uSchemeName << std::endl;
            }
        }

        std::string kSchemeName =
            convSec.lookupOrDefault<std::string>("k", "");

        if (!kSchemeName.empty())
        {
            schemes.kScheme =
                createConvectionScheme(kSchemeName);

            if (debug_)
            {
                std::cout
                    << "k convection scheme: "
                    << kSchemeName << std::endl;
            }
        }

        std::string omegaSchemeName =
            convSec.lookupOrDefault<std::string>("omega", "");

        if (!omegaSchemeName.empty())
        {
            schemes.omegaScheme =
                createConvectionScheme(omegaSchemeName);

            if (debug_)
            {
                std::cout
                    << "omega convection scheme: "
                    << omegaSchemeName << std::endl;
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
