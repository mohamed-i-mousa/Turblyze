/******************************************************************************
 * @file CFDApplication.cpp
 * @brief Top-level application driver for the CFD solver
 *****************************************************************************/

#include "CFDApplication.hpp"

#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

#include "MeshReader.hpp"
#include "MeshChecker.hpp"
#include "VtkWriter.hpp"
#include "LinearSolvers.hpp"
#include "Constraint.hpp"


// ******************************* Constructor *******************************

CFDApplication::CFDApplication(const std::string& caseFilePath)
    : caseFilePath_(caseFilePath) {}


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
    auto mesh = caseReader_->section("mesh");
    std::string meshFilePath = mesh.lookup<std::string>("file");
    checkQuality_ = mesh.lookupOrDefault<bool>("checkQuality", true);

    // Extract physical properties
    auto physicalProperties = caseReader_->section("physicalProperties");
    rho_ = physicalProperties.lookup<Scalar>("rho");
    mu_ = physicalProperties.lookup<Scalar>("mu");

    // Extract initial conditions
    auto initialConditions = caseReader_->section("initialConditions");
    initialVelocity_ = initialConditions.lookup<Vector>("U");
    initialPressure_ = initialConditions.lookup<Scalar>("p");

    // Parse convection schemes
    convectionSchemes_ = parseConvectionSchemes();

    // Extract SIMPLE parameters
    auto simple = caseReader_->section("SIMPLE");

    maxIterations_ = simple.lookup<int>("numIterations");

    convergenceTolerance_ = simple.lookup<Scalar>("convergenceTolerance");

    auto relaxFactors = simple.section("relaxationFactors");
    alphaU_ = relaxFactors.lookup<Scalar>("U");
    alphaP_ = relaxFactors.lookup<Scalar>("p");
    alphaK_ = relaxFactors.lookupOrDefault<Scalar>("k", 0.5);
    alphaOmega_ = relaxFactors.lookupOrDefault<Scalar>("omega", 0.5);

    // Extract turbulence parameters
    auto turbulence = caseReader_->section("turbulence");
    turbulenceEnabled_ = turbulence.lookup<bool>("enabled");
    turbulenceModel_ = turbulence.lookup<std::string>("model");

    // Compute turbulence initial conditions
    Scalar turbIntensity =
        turbulence.lookupOrDefault<Scalar>("turbulenceIntensity", S(0.05));
    Scalar hydraulicDiam =
        turbulence.lookupOrDefault<Scalar>("hydraulicDiameter", S(0.01));

    Scalar lTurb = S(0.07) * hydraulicDiam;
    Scalar Cmu = S(0.09);
    Scalar UMag = initialVelocity_.magnitude();

    defaultK_ = 
        std::max
        (
            S(1.5) * (turbIntensity * UMag) * (turbIntensity * UMag), S(1e-8)
        );

    defaultOmega_ = 
        std::max
        (
            std::sqrt(defaultK_) / (std::pow(Cmu, S(0.25)) * lTurb), S(1e-4)
        );

    initialK_ = initialConditions.lookupOrDefault<Scalar>("k", defaultK_);

    initialOmega_ = 
        initialConditions.lookupOrDefault<Scalar>("omega", defaultOmega_);

    // Extract output configuration
    auto outputDict = caseReader_->section("output");
    vtkOutputFilename_ = outputDict.lookup<std::string>("filename");

    // Extract constraints (optional)
    if (caseReader_->hasSection("constraints"))
    {
        auto constraintsDict = caseReader_->section("constraints");

        if (constraintsDict.hasSection("velocity"))
        {
            auto velConstraint = constraintsDict.section("velocity");

            velocityConstraintEnabled_ =
                velConstraint.lookup<bool>("enabled");

            maxVelocityConstraint_ =
                velConstraint.lookup<Scalar>("maxVelocity");
        }

        if (constraintsDict.hasSection("pressure"))
        {
            auto presConstraint = constraintsDict.section("pressure");

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

    auto mesh = caseReader_->section("mesh");
    std::string meshFilePath = mesh.lookup<std::string>("file");

    MeshReader meshReader(meshFilePath);

    nodes_ = meshReader.moveNodes();
    faces_ = meshReader.moveFaces();
    cells_ = meshReader.moveCells();
    patches_ = meshReader.moveBoundaryPatches();

    std::cout
        << "Mesh Loaded: " << nodes_.size() << " nodes, "
        << faces_.size() << " faces, " << cells_.size()
        << " cells." << std::endl;

    for (auto& face : faces_)
    {
        face.calculateGeometricProperties(nodes_);
    }
    std::cout
        << "Geometric properties calculated for faces."
        << std::endl;

    {
        std::vector<Vector> approxCentroids(cells_.size(), Vector(0,0,0));
        std::vector<std::set<size_t>> cellNodes(cells_.size());

        // Collect unique node indices for each cell
        for (const auto& face : faces_)
        {
            size_t owner = face.ownerCell();
            for (size_t f : face.nodeIndices())
            {
                cellNodes[owner].insert(f);
            }

            if (!face.isBoundary())
            {
                size_t neighbor = face.neighborCell().value();

                for (size_t n : face.nodeIndices())
                {
                    cellNodes[neighbor].insert(n);
                }
            }
        }

        // Compute centroid as average of node positions
        for (size_t i = 0; i < cells_.size(); ++i)
        {
            if (!cellNodes[i].empty())
            {
                for (size_t n : cellNodes[i])
                {
                    approxCentroids[i] += nodes_[n];
                }
                approxCentroids[i] /= S(cellNodes[i].size());
            }
        }

        // Check and correct any inverted face normals
        int flippedCount = 0;
        for (auto& face : faces_)
        {
            if (!face.isBoundary())
            {
                const Vector& ownerCell = approxCentroids[face.ownerCell()];

                const Vector& neighborCell = 
                    approxCentroids[face.neighborCell().value()];

                Vector d_PN = neighborCell - ownerCell;

                if (dot(d_PN, face.normal()) < 0)
                {
                    face.flipNormal();
                    flippedCount++;
                }
            }
        }

        if (flippedCount > 0)
        {
            std::cout
                << "Corrected " << flippedCount
                << " inverted face normals." << std::endl;
        }
    }

    for (auto& cell : cells_)
    {
        cell.calculateGeometricProperties(faces_);
    }

    std::cout
        << "Geometric properties calculated for cells."
        << std::endl;

    // Write cell geometry data for verification
    VtkWriter::writeCellGeometryData
    (
        "../outputFiles.nosync/cell_geometry_mycfdcode.txt",
        cells_
    );

    for (auto& face : faces_)
    {
        face.calculateDistanceProperties(cells_);
    }
    std::cout
        << "Distance properties calculated for faces."
        << std::endl;

    // Check mesh quality if requested
    if (checkQuality_)
    {
        MeshChecker meshChecker(nodes_, faces_, cells_);
        meshChecker.check();
    }
}


// ************************* setupBoundaryConditions *************************

void CFDApplication::setupBoundaryConditions()
{
    std::cout
        << std::endl << "--- 2. Setting Boundary Conditions ---" << std::endl;

    for (const auto& patch : patches_)
    {
        bcManager_.addPatch(patch);
    }

    // Load boundary conditions from case file
    auto BCs = caseReader_->section("boundaryConditions");

    std::set<std::string> omegaWallFunctionPatches;
    std::set<std::string> legacyKFixedZeroPatches;

    // Process velocity boundary conditions
    if (BCs.hasSection("U"))
    {
        auto velocityBCs = BCs.section("U");

        for (const auto& patchName : velocityBCs.sectionNames())
        {
            auto patchBC = velocityBCs.section(patchName);
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                Vector value = patchBC.lookup<Vector>("value");
                bcManager_.setFixedValue(patchName, "U", value);
                
                bcManager_.setFixedValue(patchName, "Ux", value.x());
                bcManager_.setFixedValue(patchName, "Uy", value.y());
                bcManager_.setFixedValue(patchName, "Uz", value.z());
            }
            else if (bcType == "noSlip")
            {
                bcManager_.setNoSlip(patchName, "U");

                bcManager_.setFixedValue(patchName, "Ux", Scalar(0));
                bcManager_.setFixedValue(patchName, "Uy", Scalar(0));
                bcManager_.setFixedValue(patchName, "Uz", Scalar(0));
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "U");

                bcManager_.setZeroGradient(patchName, "Ux");
                bcManager_.setZeroGradient(patchName, "Uy");
                bcManager_.setZeroGradient(patchName, "Uz");
            }
        }
    }

    // Process pressure boundary conditions
    if (BCs.hasSection("p"))
    {
        auto pressureBCs = BCs.section("p");

        for (const auto& patchName : pressureBCs.sectionNames())
        {
            auto patchBC = pressureBCs.section(patchName);
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                Scalar value = patchBC.lookup<Scalar>("value");
                bcManager_.setFixedValue(patchName, "p", value);
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "p");
            }
        }

        // Derive pressure correction BCs from pressure BCs
        for (const auto& patchName
            : pressureBCs.sectionNames())
        {
            auto patchBC = pressureBCs.section(patchName);
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

    // Process wall distance field boundary conditions
    if (BCs.hasSection("phi_wall"))
    {
        auto phiWallBCs = BCs.section("phi_wall");

        for (const auto& patchName : phiWallBCs.sectionNames())
        {
            auto patchBC = phiWallBCs.section(patchName);
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                Scalar value = patchBC.lookup<Scalar>("value");
                bcManager_.setFixedValue(patchName, "phi_wall", value);
            }
        }
    }

    // Process turbulent kinetic energy BCs
    if (BCs.hasSection("k"))
    {
        auto kBCs = BCs.section("k");

        for (const auto& patchName : kBCs.sectionNames())
        {
            auto patchBC = kBCs.section(patchName);
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                std::string valStr = patchBC.lookup<std::string>("value");
                Scalar value = 
                    (valStr == "calculated") ? 
                        defaultK_
                      : patchBC.lookup<Scalar>("value");

                bcManager_.setFixedValue(patchName, "k", value);

                if (std::abs(value) <= S(1e-14))
                {
                    legacyKFixedZeroPatches.insert(patchName);
                }
            }
            else if (bcType == "kWallFunction")
            {
                bcManager_.setKWallFunction(patchName, "k");
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "k");
            }
        }
    }

    // Process specific dissipation rate BCs
    if (BCs.hasSection("omega"))
    {
        auto omegaBCs = BCs.section("omega");

        for (const auto& patchName : omegaBCs.sectionNames())
        {
            auto patchBC = omegaBCs.section(patchName);
            std::string bcType = patchBC.lookup<std::string>("type");

            if (bcType == "fixedValue")
            {
                std::string valStr = patchBC.lookup<std::string>("value");
                Scalar value =
                    (valStr == "calculated") ?
                    defaultOmega_
                  : patchBC.lookup<Scalar>("value");
                
                  bcManager_.setFixedValue(patchName, "omega", value);
            }
            else if (bcType == "omegaWallFunction")
            {
                bcManager_.setOmegaWallFunction(patchName, "omega");
                omegaWallFunctionPatches.insert(patchName);
            }
            else if (bcType == "zeroGradient")
            {
                bcManager_.setZeroGradient(patchName, "omega");
            }
        }
    }

    // Optional turbulent viscosity BC section
    if (BCs.hasSection("nut"))
    {
        auto nutBCs = BCs.section("nut");

        for (const auto& patchName : nutBCs.sectionNames())
        {
            auto patchBC = nutBCs.section(patchName);
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
            else if (bcType == "nutkWallFunction")
            {
                bcManager_.setNutWallFunction(patchName, "nut");
            }
        }
    }

    bcManager_.printSummary();
}


// ***************************** configureSolver *****************************

void CFDApplication::configureSolver()
{
    std::cout
        << std::endl << "--- 3. Initializing SIMPLE Solver ---" << std::endl;

    gradScheme_ = 
        std::make_unique<GradientScheme>(faces_, cells_, bcManager_);

    solver_ = 
        std::make_unique<SIMPLE>
        (
            faces_,
            cells_,
            bcManager_,
            *gradScheme_,
            convectionSchemes_
        );

    // Configure turbulence modeling
    solver_->enableTurbulenceModeling(turbulenceEnabled_);
    if (turbulenceEnabled_)
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
            initialOmega_
        );

    solver_->setPhysicalProperties(rho_, mu_);

    // Configure SIMPLE parameters
    solver_->setRelaxationFactors(alphaU_, alphaP_, alphaK_, alphaOmega_);
    solver_->setConvergenceTolerance(convergenceTolerance_);
    solver_->setMaxIterations(maxIterations_);

    // Configure linear solvers from case file
    if (caseReader_->hasSection("linearSolvers"))
    {
        auto solvers = caseReader_->section("linearSolvers");

        // Momentum solver (U)
        if (solvers.hasSection("U"))
        {
            auto USection = solvers.section("U");
            LinearSolver momentumSolver
            (
                "momentum",
                USection.lookupOrDefault<Scalar>("tolerance", S(1e-8)),
                USection.lookupOrDefault<int>("maxIter", 1000)
            );
            momentumSolver.setComputeResiduals
            (
                USection.lookupOrDefault<bool>("computeResiduals", false)
            );
            solver_->setMomentumSolver(momentumSolver);

            std::cout
                << "  Momentum solver: tolerance="
                << USection.lookupOrDefault<Scalar>("tolerance", S(1e-8))
                << ", maxIter="
                << USection.lookupOrDefault<int>("maxIter", 1000)
                << std::endl;
        }

        // Pressure solver (p)
        if (solvers.hasSection("p"))
        {
            auto pSection = solvers.section("p");
            LinearSolver pressureSolver
            (
                "pCorr",
                pSection.lookupOrDefault<Scalar>("tolerance", S(1e-6)),
                pSection.lookupOrDefault<int>("maxIter", 1000),
                pSection.lookupOrDefault<Scalar>("relTol", S(0.05))
            );

            Scalar initialShift =
                pSection.lookupOrDefault<Scalar>("initialShift", S(1e-2));

            pressureSolver.setICParameters(initialShift);

            pressureSolver.setComputeResiduals
            (
                pSection.lookupOrDefault<bool>("computeResiduals", false)
            );

            solver_->setPressureSolver(pressureSolver);

            std::cout
                << "  Pressure solver: tolerance="
                << pSection.lookupOrDefault<Scalar>("tolerance", S(1e-6))
                << ", relTol="
                << pSection.lookupOrDefault<Scalar>(
                    "relTol", S(0.05))
                << ", maxIter="
                << pSection.lookupOrDefault<int>(
                    "maxIter", 1000)
                << std::endl;
        }

        // Turbulence solvers (k, omega)
        if (turbulenceEnabled_)
        {
            LinearSolver kSolver
            (
                "k",
                solvers.hasSection("k") ? 
                solvers.section("k").lookupOrDefault<Scalar>
                (
                    "tolerance",
                    S(1e-8)
                )
              : S(1e-8),

                solvers.hasSection("k") ? 
                solvers.section("k").lookupOrDefault<int>("maxIter", 1000)
              : 1000,

                solvers.hasSection("k") ?
                solvers.section("k").lookupOrDefault<Scalar>("relTol", S(0.1))
              : S(0.1)
            );
            if (solvers.hasSection("k"))
            {
                kSolver.setComputeResiduals
                (
                    solvers.section("k").lookupOrDefault<bool>
                    (
                        "computeResiduals",
                        false
                    )
                );
            }

            LinearSolver omegaSolver
            (
                "omega",
                solvers.hasSection("omega") ?
                solvers.section("omega").lookupOrDefault<Scalar>
                (
                    "tolerance", 
                    S(1e-8)
                )
              : S(1e-8),

                solvers.hasSection("omega") ?
                solvers.section("omega").lookupOrDefault<int>("maxIter", 1000)
              : 1000,

                solvers.hasSection("omega") ?
                solvers.section("omega").lookupOrDefault<Scalar>
                (
                    "relTol", 
                    S(0.1)
                )
              : S(0.1)
            );
            if (solvers.hasSection("omega"))
            {
                omegaSolver.setComputeResiduals
                (
                    solvers.section("omega").lookupOrDefault<bool>
                    (
                        "computeResiduals",
                        false
                    )
                );
            }

            solver_->setTurbulenceSolvers(kSolver, omegaSolver);
        }
    }

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

    // Configure field constraints
    Constraint* constraintSystem = solver_->getConstraintSystem();

    if (constraintSystem)
    {
        if (velocityConstraintEnabled_)
        {
            constraintSystem->setVelocityConstraints(maxVelocityConstraint_);
            std::cout
                << "Velocity constraint enabled: max = "
                << maxVelocityConstraint_ << " m/s"
                << std::endl;
        }
        if (pressureConstraintEnabled_)
        {
            constraintSystem->setPressureConstraints
            (
                minPressureConstraint_,
                maxPressureConstraint_
            );

            std::cout
                << "Pressure constraint enabled: ["
                << minPressureConstraint_
                << ", " << maxPressureConstraint_
                << "] Pa" << std::endl;
        }
        constraintSystem->enableConstraints
        (
            velocityConstraintEnabled_,
            pressureConstraintEnabled_
        );
    }
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

    const VectorField& velocity = solver_->getVelocity();
    const ScalarField& pressure = solver_->getPressure();

    std::cout
        << "Solution extracted." << std::endl;

    std::cout
         << std::endl << "--- 6. Post-Processing Results ---" << std::endl;

    // Calculate velocity magnitude
    ScalarField velocityMagnitude =
        VtkWriter::computeVelocityMagnitude(velocity);

    // Print statistics
    Scalar maximumVelocity = 0.0;
    Scalar averageVelocity = 0.0;
    Scalar maximumPressure = pressure[0];
    Scalar minimumPressure = pressure[0];

    for (size_t i = 0; i < velocity.size(); ++i)
    {
        Scalar vmag = velocityMagnitude[i];
        maximumVelocity = std::max(maximumVelocity, vmag);
        averageVelocity += vmag;

        maximumPressure = std::max(maximumPressure, pressure[i]);
        minimumPressure = std::min(minimumPressure, pressure[i]);
    }
    averageVelocity /= velocity.size();

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

    const VectorField& velocity = solver_->getVelocity();
    const ScalarField& pressure = solver_->getPressure();

    // Extract turbulence fields if available
    const ScalarField* kField = solver_->getTurbulentKineticEnergy();
    const ScalarField* omegaField = solver_->getSpecificDissipationRate();
    const ScalarField* nutField = solver_->getTurbulentViscosity();
    const ScalarField* wallDistField = solver_->getWallDistance();

    // Calculate velocity magnitude
    ScalarField velocityMagnitude =
        VtkWriter::computeVelocityMagnitude(velocity);

    // Prepare scalar fields for export
    std::map<std::string, const ScalarField*>
    scalarFieldsToVtk;

    scalarFieldsToVtk["pressure"] = &pressure;

    scalarFieldsToVtk["velocityMagnitude"] = &velocityMagnitude;

    // Add turbulence fields if available
    if
    (
        turbulenceEnabled_
     && kField
     && omegaField
     && nutField
     && wallDistField
    )
    {
        scalarFieldsToVtk["k"] = kField;
        scalarFieldsToVtk["omega"] = omegaField;
        scalarFieldsToVtk["nut"] = nutField;
        scalarFieldsToVtk["wallDistance"] = wallDistField;
    }

    // Prepare vector fields for export
    std::map<std::string, const VectorField*>
    vectorFieldsToVtk;
    
    vectorFieldsToVtk["velocity"] = &velocity;

    // Determine VTU filename
    std::string vtuFilename = vtkOutputFilename_;
    size_t pos = vtuFilename.rfind(".vtu");
    if (pos != std::string::npos)
    {
        vtuFilename.replace(pos, 4, ".vtu");
    }
    else
    {
        if (vtuFilename.find(".vtu") == std::string::npos)
        {
            vtuFilename += ".vtu";
        }
    }

    std::cout
        << std::endl << "Exporting results to VTK UnstructuredGrid..."
        << std::endl;

    VtkWriter::writeVtkUnstructuredGrid
    (
        vtuFilename,
        nodes_,
        cells_,
        faces_,
        scalarFieldsToVtk,
        vectorFieldsToVtk
    );

    std::cout
        << std::endl << "=== CFD Results Exported Successfully ==="
        << std::endl;
    std::cout
        << "File: " << vtuFilename << std::endl;
}


// ************************* createConvectionScheme ***********************

std::unique_ptr<ConvectionScheme>
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
        throw   std::runtime_error
                (
                    "Unknown convection scheme: " + name
                );
    }
}


// ************************ parseConvectionSchemes ************************

ConvectionSchemes CFDApplication::parseConvectionSchemes()
{
    auto schemesDict = caseReader_->section("numericalSchemes");

    ConvectionSchemes schemes;

    if (schemesDict.hasSection("convection"))
    {
        auto convSec = schemesDict.section("convection");

        std::string defaultSchemeName = 
            convSec.lookupOrDefault<std::string>("default", "Upwind");

        schemes.defaultScheme = createConvectionScheme(defaultSchemeName);

        std::cout
            << "Default convection scheme: "
            << defaultSchemeName << std::endl;

        // Per-equation overrides (fall back to default if unset)
        std::string uSchemeName =
            convSec.lookupOrDefault<std::string>("U", "");

        if (!uSchemeName.empty())
        {
            schemes.momentumScheme = createConvectionScheme(uSchemeName);
            
            std::cout
                << "Momentum convection scheme: "
                << uSchemeName << std::endl;
        }

        std::string kSchemeName =
            convSec.lookupOrDefault<std::string>("k", "");

        if (!kSchemeName.empty())
        {
            schemes.kScheme = createConvectionScheme(kSchemeName);
            
            std::cout
                << "k convection scheme: "
                << kSchemeName << std::endl;
        }

        std::string omegaSchemeName =
            convSec.lookupOrDefault<std::string>("omega", "");

        if (!omegaSchemeName.empty())
        {
            schemes.omegaScheme = createConvectionScheme(omegaSchemeName);
            
            std::cout
                << "omega convection scheme: "
                << omegaSchemeName << std::endl;
        }
    }
    else
    {
        throw   std::runtime_error
                (
                    "Missing 'convection' sub-section "
                    "in numericalSchemes"
                );
    }

    return schemes;
}