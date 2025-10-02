/******************************************************************************
 * @file main.cpp
 * @brief Main entry point for the 3D incompressible CFD solver
 * 
 * This file contains the main function that orchestrates the entire CFD
 * simulation process using the SIMPLE algorithm. It handles mesh reading, 
 * boundary condition setup, solver setup, solution computation, 
 * and results export.
 * 
 * The solver implements:
 * - 3D incompressible Navier-Stokes equations
 * - SIMPLE algorithm for pressure-velocity coupling
 * - k-omega SST turbulence modeling
 *
 * @author Mohamed Mousa
 * @date 2025
 *****************************************************************************/

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <map>
#include <algorithm>
#include <chrono>

// Core data structures
#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"

// Boundary conditions
#include "BoundaryPatch.hpp"
#include "BoundaryData.hpp"
#include "BoundaryConditions.hpp"

// Numerics and Solver
#include "GradientScheme.hpp"
#include "ConvectionScheme.hpp"
#include "Matrix.hpp"
#include "LinearSolvers.hpp"
#include "SIMPLE.hpp"
#include "Constraint.hpp"

// I/O
#include "MeshReader.hpp"
#include "VtkWriter.hpp"
#include "checkMesh.hpp"

// Setup
#include "SetupReader.hpp"

int main(int argc, char* argv[])
{
    // Start timing the total execution
    auto start_time = std::chrono::high_resolution_clock::now();

    std::cout   << "--- Welcome to the CFD Solver with SIMPLE Algorithm ---"
                << std::endl;

    std::cout   << "Running with precision: " << SCALAR_MODE
                << std::endl;

    std::cout   << std::fixed << std::setprecision(6);

    try
    {
        /**********************************************************************
         * ----------------------- 0. LOAD SETUP -----------------------------
         *********************************************************************/

        std::string setupFile = "../defaultSetup";

        if (argc > 1)
        {
            setupFile = argv[1];

            std::cout   << "Using setup file: " << setupFile
                        << std::endl;
        }
        else
        {
            std::cout   << "Using default setup: " << setupFile
                        << std::endl;
        }

        SetupReader setup(setupFile);

        // Extract mesh setup
        auto mesh = setup.section("mesh");
        std::string meshFilePath = mesh.lookup<std::string>("file");
        bool checkQuality = mesh.lookupOrDefault<bool>("checkQuality", true);

        // Extract physical properties
        auto physProps = setup.section("physicalProperties");
        const Scalar rho = physProps.lookup<Scalar>("rho");
        const Scalar mu = physProps.lookup<Scalar>("mu");

        // Extract initial conditions
        auto initConds = setup.section("initialConditions");
        Vector initialVelocity = initConds.lookup<Vector>("U");
        Scalar initialPressure = initConds.lookup<Scalar>("p");

        // Extract numerical schemes
        auto schemes = setup.section("numericalSchemes");
        std::string convectionScheme = schemes.lookup<std::string>("convection");

        // Extract SIMPLE parameters
        auto simpleDict = setup.section("SIMPLE");
        int maxIterations = simpleDict.lookup<int>("nIterations");
        Scalar convergenceTolerance = simpleDict.lookup<Scalar>("convergenceTolerance");

        auto relaxFactors = simpleDict.section("relaxationFactors");
        Scalar alphaU = relaxFactors.lookup<Scalar>("U");
        Scalar alphaP = relaxFactors.lookup<Scalar>("p");
        Scalar alphaK = relaxFactors.lookupOrDefault<Scalar>("k", 0.5);
        Scalar alphaOmega = relaxFactors.lookupOrDefault<Scalar>("omega", 0.5);

        // Extract turbulence parameters
        auto turbDict = setup.section("turbulence");
        bool turbulenceEnabled = turbDict.lookup<bool>("enabled");
        std::string turbulenceModel = turbDict.lookup<std::string>("model");

        // Extract output setup
        auto outputDict = setup.section("output");
        std::string vtkOutputFilename = outputDict.lookup<std::string>("filename");

        // Extract constraints (optional)
        bool velocityConstraintEnabled = false;
        Scalar maxVelocity = 0.3;
        bool pressureConstraintEnabled = false;
        Scalar minPressure = -0.05;
        Scalar maxPressure = 0.05;

        if (setup.hasSection("constraints")) {
            auto constraintsDict = setup.section("constraints");
            if (constraintsDict.hasSection("velocity")) {
                auto velConstraint = constraintsDict.section("velocity");
                velocityConstraintEnabled = velConstraint.lookup<bool>("enabled");
                maxVelocity = velConstraint.lookup<Scalar>("maxVelocity");
            }
            if (constraintsDict.hasSection("pressure")) {
                auto presConstraint = constraintsDict.section("pressure");
                pressureConstraintEnabled = presConstraint.lookup<bool>("enabled");
                minPressure = presConstraint.lookup<Scalar>("minPressure");
                maxPressure = presConstraint.lookup<Scalar>("maxPressure");
            }
        }

        /**********************************************************************
         * -------------------------- 1. MESH SETUP ---------------------------
         *********************************************************************/

        std::cout << "\n--- 1. Reading and Preparing Mesh ---" << std::endl;

        std::vector<Vector> allNodes;
        std::vector<Face> allFaces;
        std::vector<Cell> allCells;
        std::vector<BoundaryPatch> allBoundaryPatches;

        readMshFile
        (
            meshFilePath,
            allNodes,
            allFaces,
            allCells,
            allBoundaryPatches
        );

        std::cout << "Mesh Loaded: " << allNodes.size() << " nodes, "
                  << allFaces.size() << " faces, " << allCells.size()
                  << " cells." << std::endl;

        for (auto& face : allFaces) {
            face.calculateGeometricProperties(allNodes);
        }
        std::cout << "Geometric properties calculated for faces." << std::endl;

        for (auto& cell : allCells) {
            cell.calculateGeometricProperties(allFaces);
        }
        std::cout << "Geometric properties calculated for cells." << std::endl;

        for (auto& face : allFaces) {
            face.calculateDistanceProperties(allCells);
        }
        std::cout << "Distance properties calculated for faces." << std::endl;

        // Check mesh quality if requested
        if (checkQuality) {
            checkMesh(allFaces, allCells);
        }

        /**********************************************************************
         * ---------------------- 2. BOUNDARY CONDITIONS ----------------------
         *********************************************************************/

        std::cout << "\n--- 2. Setting up Boundary Conditions ---" << std::endl;

        BoundaryConditions bcManager;
        for (const auto& patch : allBoundaryPatches) {
            bcManager.addPatch(patch);
        }

        // Load boundary conditions from setup
        auto bcDict = setup.section("boundaryConditions");

        // Process velocity boundary conditions
        if (bcDict.hasSection("U")) {
            auto velocityBCs = bcDict.section("U");

            for (const auto& patchName : velocityBCs.sectionNames()) {
                auto patchBC = velocityBCs.section(patchName);
                std::string bcType = patchBC.lookup<std::string>("type");

                if (bcType == "fixedValue") {
                    Vector value = patchBC.lookup<Vector>("value");
                    bcManager.setFixedValue(patchName, "U", value);
                } else if (bcType == "noSlip") {
                    bcManager.setNoSlip(patchName, "U");
                } else if (bcType == "zeroGradient") {
                    bcManager.setZeroGradient(patchName, "U");
                }
            }
        }

        // Process pressure boundary conditions
        if (bcDict.hasSection("p")) {
            auto pressureBCs = bcDict.section("p");

            for (const auto& patchName : pressureBCs.sectionNames()) {
                auto patchBC = pressureBCs.section(patchName);
                std::string bcType = patchBC.lookup<std::string>("type");

                if (bcType == "fixedValue") {
                    Scalar value = patchBC.lookup<Scalar>("value");
                    bcManager.setFixedValue(patchName, "p", value);
                } else if (bcType == "zeroGradient") {
                    bcManager.setZeroGradient(patchName, "p");
                }
            }
        }

        // Process wall distance field boundary conditions (for turbulence)
        if (bcDict.hasSection("phi_wall")) {
            auto phiWallBCs = bcDict.section("phi_wall");

            for (const auto& patchName : phiWallBCs.sectionNames()) {
                auto patchBC = phiWallBCs.section(patchName);
                std::string bcType = patchBC.lookup<std::string>("type");

                if (bcType == "fixedValue") {
                    Scalar value = patchBC.lookup<Scalar>("value");
                    bcManager.setFixedValue(patchName, "phi_wall", value);
                }
            }
        }

        // Process turbulent kinetic energy boundary conditions
        if (bcDict.hasSection("k")) {
            auto kBCs = bcDict.section("k");

            for (const auto& patchName : kBCs.sectionNames()) {
                auto patchBC = kBCs.section(patchName);
                std::string bcType = patchBC.lookup<std::string>("type");

                if (bcType == "fixedValue") {
                    Scalar value = patchBC.lookup<Scalar>("value");
                    bcManager.setFixedValue(patchName, "k", value);
                } else if (bcType == "zeroGradient") {
                    bcManager.setZeroGradient(patchName, "k");
                }
            }
        }

        // Process specific dissipation rate boundary conditions
        if (bcDict.hasSection("omega")) {
            auto omegaBCs = bcDict.section("omega");

            for (const auto& patchName : omegaBCs.sectionNames()) {
                auto patchBC = omegaBCs.section(patchName);
                std::string bcType = patchBC.lookup<std::string>("type");

                if (bcType == "fixedValue") {
                    Scalar value = patchBC.lookup<Scalar>("value");
                    bcManager.setFixedValue(patchName, "omega", value);
                } else if (bcType == "zeroGradient") {
                    bcManager.setZeroGradient(patchName, "omega");
                }
            }
        }

        bcManager.printSummary();

        /**********************************************************************
         * --------------------- 3. SIMPLE SOLVER SETUP -----------------------
         *********************************************************************/

        std::cout << "\n--- 3. Initializing SIMPLE Solver ---" << std::endl;

        // Select convection scheme based on setup
        ConvectionScheme* selectedScheme = nullptr;
        std::unique_ptr<UpwindScheme> uds;
        std::unique_ptr<CentralDifferenceScheme> cds;
        std::unique_ptr<SecondOrderUpwindScheme> sous;

        if (convectionScheme == "Upwind") {
            uds = std::make_unique<UpwindScheme>();
            selectedScheme = uds.get();
            std::cout << "Using Upwind convection scheme" << std::endl;
        } else if (convectionScheme == "CentralDifference") {
            cds = std::make_unique<CentralDifferenceScheme>();
            selectedScheme = cds.get();
            std::cout << "Using Central Difference convection scheme" << std::endl;
        } else if (convectionScheme == "SecondOrderUpwind") {
            sous = std::make_unique<SecondOrderUpwindScheme>();
            selectedScheme = sous.get();
            std::cout << "Using Second Order Upwind convection scheme" << std::endl;
        } else {
            throw std::runtime_error("Unknown convection scheme: " + convectionScheme);
        }

        GradientScheme gradScheme;

        SIMPLE simpleSolver(allFaces, allCells, bcManager, gradScheme, *selectedScheme);

        // Setup turbulence modeling BEFORE initialization
        // (turbulence model is created inside initialize() if enabled)
        simpleSolver.enableTurbulenceModeling(turbulenceEnabled);
        if (turbulenceEnabled) {
            std::cout << "Turbulence modeling enabled: " << turbulenceModel << std::endl;
        }

        // Initialize velocity and pressure fields
        // (this creates the turbulence model object if turbulence is enabled)
        simpleSolver.initialize(initialVelocity, initialPressure);
        simpleSolver.setPhysicalProperties(rho, mu);

        // Setup SIMPLE parameters
        simpleSolver.setRelaxationFactors(alphaU, alphaP, alphaK, alphaOmega);
        simpleSolver.setConvergenceTolerance(convergenceTolerance);
        simpleSolver.setMaxIterations(maxIterations);

        std::cout << "SIMPLE parameters:" << std::endl;
        std::cout << "  Max iterations: " << maxIterations << std::endl;
        std::cout << "  Convergence tolerance: " << convergenceTolerance << std::endl;
        std::cout << "  Velocity relaxation: " << alphaU << std::endl;
        std::cout << "  Pressure relaxation: " << alphaP << std::endl;
        if (turbulenceEnabled) {
            std::cout << "  k relaxation: " << alphaK << std::endl;
            std::cout << "  omega relaxation: " << alphaOmega << std::endl;
        }

        // Setup field constraints
        Constraint* constraintSystem = simpleSolver.getConstraintSystem();
        if (constraintSystem) {
            if (velocityConstraintEnabled) {
                constraintSystem->setVelocityConstraints(maxVelocity);
                std::cout << "Velocity constraint enabled: max = " << maxVelocity << " m/s" << std::endl;
            }
            if (pressureConstraintEnabled) {
                constraintSystem->setPressureConstraints(minPressure, maxPressure);
                std::cout << "Pressure constraint enabled: [" << minPressure << ", " << maxPressure << "] Pa" << std::endl;
            }
            constraintSystem->enableConstraints(velocityConstraintEnabled, pressureConstraintEnabled);
        }

        /**********************************************************************
         * ------------------- 4. SOLVE STEADY-STATE FLOW ---------------------
         *********************************************************************/

        std::cout << "\n--- 4. Solving Steady-State Flow with SIMPLE ---" << std::endl;

        simpleSolver.solve();

        /**********************************************************************
         * ------------------- 5. EXTRACT SOLUTION FIELDS ---------------------
         *********************************************************************/

        std::cout << "\n--- 5. Extracting Solution Fields ---" << std::endl;

        const VectorField& velocity = simpleSolver.getVelocity();
        const ScalarField& pressure = simpleSolver.getPressure();

        // Extract turbulence fields if available
        const ScalarField* k_field = simpleSolver.getTurbulentKineticEnergy();
        const ScalarField* omega_field = simpleSolver.getSpecificDissipationRate();
        const ScalarField* nu_t_field = simpleSolver.getTurbulentViscosity();
        const ScalarField* wallDist_field = simpleSolver.getWallDistance();

        // Extract 3D velocity components
        ScalarField U_x = simpleSolver.getVelocityX();
        ScalarField U_y = simpleSolver.getVelocityY();
        ScalarField U_z = simpleSolver.getVelocityZ();

        std::cout << "Solution extracted." << std::endl;

        /**********************************************************************
         * ----------------------- 6. POST-PROCESSING -------------------------
         *********************************************************************/

        std::cout << "\n--- 6. Post-Processing Results ---" << std::endl;

        // Calculate velocity magnitude
        ScalarField velocityMagnitude("velocityMagnitude", velocity.size());

        for (size_t i = 0; i < velocity.size(); ++i) {
            velocityMagnitude[i] = velocity[i].magnitude();
        }

        // Print statistics
        Scalar maxVel = 0.0, avgVel = 0.0;
        Scalar maxPress = pressure[0], minPress = pressure[0];

        for (size_t i = 0; i < velocity.size(); ++i) {
            Scalar vmag = velocityMagnitude[i];
            maxVel = std::max(maxVel, vmag);
            avgVel += vmag;

            maxPress = std::max(maxPress, pressure[i]);
            minPress = std::min(minPress, pressure[i]);
        }
        avgVel /= velocity.size();

        std::cout << "Flow Statistics:" << std::endl;
        std::cout << "  Max velocity magnitude: " << maxVel << " m/s" << std::endl;
        std::cout << "  Average velocity magnitude: " << avgVel << " m/s" << std::endl;
        std::cout << "  Pressure range: [" << minPress << ", " << maxPress << "] Pa" << std::endl;

        /**********************************************************************
         * ------------------------ 7. EXPORT RESULTS -------------------------
         *********************************************************************/

        std::cout << "\n--- 7. Exporting Results to VTK ---" << std::endl;

        // Prepare scalar fields for export
        std::map<std::string, const ScalarField*> scalarFieldsToVtk;
        scalarFieldsToVtk["pressure"] = &pressure;
        scalarFieldsToVtk["velocityMagnitude"] = &velocityMagnitude;
        scalarFieldsToVtk["U_x"] = &U_x;
        scalarFieldsToVtk["U_y"] = &U_y;
        scalarFieldsToVtk["U_z"] = &U_z;

        // Add turbulence fields if available
        if (turbulenceEnabled && k_field && omega_field && nu_t_field && wallDist_field) {
            scalarFieldsToVtk["k"] = k_field;
            scalarFieldsToVtk["omega"] = omega_field;
            scalarFieldsToVtk["nu_t"] = nu_t_field;
            scalarFieldsToVtk["wallDistance"] = wallDist_field;
        }

        // Prepare vector fields for export
        std::map<std::string, const VectorField*> vectorFieldsToVtk;
        vectorFieldsToVtk["velocity"] = &velocity;

        // Export VTK UnstructuredGrid file (VTU) - the primary output format
        std::string vtuFilename = vtkOutputFilename;
        size_t pos = vtuFilename.rfind(".vtp");
        if (pos != std::string::npos)
        {
            vtuFilename.replace(pos, 4, ".vtu");
        }
        else
        {
            // Check if already has .vtu extension
            if (vtuFilename.find(".vtu") == std::string::npos)
            {
                vtuFilename += ".vtu";
            }
        }

        std::cout   << "\nExporting results to VTK UnstructuredGrid..."
                    << std::endl;
        VtkWriter::writeVtkUnstructuredGrid(
            vtuFilename,
            allNodes,
            allCells,
            allFaces,
            scalarFieldsToVtk,
            vectorFieldsToVtk
        );

        std::cout << "\n=== CFD Results Exported Successfully ===" << std::endl;
        std::cout << "File: " << vtuFilename << std::endl;

        /**********************************************************************
         * -------------------------- 8. FINALIZE ------------------------------
         *********************************************************************/

        auto end_time = std::chrono::high_resolution_clock::now();

        auto duration = 
            std::chrono::duration_cast<std::chrono::seconds>
            (end_time - start_time);

        std::cout << "\n--- Simulation Complete ---" << std::endl;
        std::cout << "Total execution time: " << duration.count() << " seconds" << std::endl;
        std::cout << "Setup file used: " << setupFile << std::endl;
    }
    catch (const std::exception& e) 
    {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}