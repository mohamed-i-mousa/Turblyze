/******************************************************************************
 * @file main.cpp
 * @brief Main entry point for the 3D incompressible CFD solver with config file support
 *
 * This version uses DictionaryReader to load configuration from a dictionary file
 * instead of hard-coding parameters.
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

// Configuration
#include "DictionaryReader.hpp"

int main(int argc, char* argv[])
{
    // Start timing the total execution
    auto start_time = std::chrono::high_resolution_clock::now();

    std::cout << "--- Welcome to the CFD Solver with SIMPLE Algorithm ---" << std::endl;
    std::cout << "Running with precision: " << SCALAR_MODE << std::endl;
    std::cout << std::fixed << std::setprecision(6);

    try
    {
        /**********************************************************************
         * -------------------- 0. LOAD CONFIGURATION ------------------------
         *********************************************************************/

        std::string configFile = "inputSettings";  // Default configuration file
        if (argc > 1) {
            configFile = argv[1];
            std::cout << "Using configuration file: " << configFile << std::endl;
        } else {
            std::cout << "Using default configuration: " << configFile << std::endl;
        }

        DictionaryReader config(configFile);

        // Extract mesh configuration
        auto meshDict = config.subDict("mesh");
        std::string meshFilePath = meshDict.lookup<std::string>("file");
        bool checkQuality = meshDict.lookupOrDefault<bool>("checkQuality", true);

        // Extract physical properties
        auto physProps = config.subDict("physicalProperties");
        const Scalar rho = physProps.lookup<Scalar>("rho");
        const Scalar mu = physProps.lookup<Scalar>("mu");

        // Extract initial conditions
        auto initConds = config.subDict("initialConditions");
        Vector initialVelocity = initConds.lookup<Vector>("U");
        Scalar initialPressure = initConds.lookup<Scalar>("p");

        // Extract numerical schemes
        auto schemes = config.subDict("numericalSchemes");
        std::string convectionScheme = schemes.lookup<std::string>("convection");

        // Extract SIMPLE settings
        auto simpleDict = config.subDict("SIMPLE");
        int maxIterations = simpleDict.lookup<int>("nIterations");
        Scalar convergenceTolerance = simpleDict.lookup<Scalar>("convergenceTolerance");

        auto relaxFactors = simpleDict.subDict("relaxationFactors");
        Scalar alphaU = relaxFactors.lookup<Scalar>("U");
        Scalar alphaP = relaxFactors.lookup<Scalar>("p");

        // Extract turbulence settings
        auto turbDict = config.subDict("turbulence");
        bool turbulenceEnabled = turbDict.lookup<bool>("enabled");
        std::string turbulenceModel = turbDict.lookup<std::string>("model");

        // Extract output configuration
        auto outputDict = config.subDict("output");
        std::string vtkOutputFilename = outputDict.lookup<std::string>("filename");

        // Extract constraints (optional)
        bool velocityConstraintEnabled = false;
        Scalar maxVelocity = 0.3;
        bool pressureConstraintEnabled = false;
        Scalar minPressure = -0.05;
        Scalar maxPressure = 0.05;

        if (config.foundSubDict("constraints")) {
            auto constraintsDict = config.subDict("constraints");
            if (constraintsDict.foundSubDict("velocity")) {
                auto velConstraint = constraintsDict.subDict("velocity");
                velocityConstraintEnabled = velConstraint.lookup<bool>("enabled");
                maxVelocity = velConstraint.lookup<Scalar>("maxVelocity");
            }
            if (constraintsDict.foundSubDict("pressure")) {
                auto presConstraint = constraintsDict.subDict("pressure");
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

        readMshFile(
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

        // Load boundary conditions from configuration
        auto bcDict = config.subDict("boundaryConditions");

        // Process velocity boundary conditions
        if (bcDict.foundSubDict("U")) {
            auto velocityBCs = bcDict.subDict("U");

            for (const auto& patchName : velocityBCs.subDictNames()) {
                auto patchBC = velocityBCs.subDict(patchName);
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
        if (bcDict.foundSubDict("p")) {
            auto pressureBCs = bcDict.subDict("p");

            for (const auto& patchName : pressureBCs.subDictNames()) {
                auto patchBC = pressureBCs.subDict(patchName);
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
        if (bcDict.foundSubDict("phi_wall")) {
            auto phiWallBCs = bcDict.subDict("phi_wall");

            for (const auto& patchName : phiWallBCs.subDictNames()) {
                auto patchBC = phiWallBCs.subDict(patchName);
                std::string bcType = patchBC.lookup<std::string>("type");

                if (bcType == "fixedValue") {
                    Scalar value = patchBC.lookup<Scalar>("value");
                    bcManager.setFixedValue(patchName, "phi_wall", value);
                }
            }
        }

        bcManager.printSummary();

        /**********************************************************************
         * --------------------- 3. SIMPLE SOLVER SETUP -----------------------
         *********************************************************************/

        std::cout << "\n--- 3. Initializing SIMPLE Solver ---" << std::endl;

        // Select convection scheme based on configuration
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

        // Configure turbulence modeling
        simpleSolver.enableTurbulenceModeling(turbulenceEnabled);
        if (turbulenceEnabled) {
            std::cout << "Turbulence modeling enabled: " << turbulenceModel << std::endl;
        }

        // Initialize velocity and pressure fields
        simpleSolver.initialize(initialVelocity, initialPressure);
        simpleSolver.setPhysicalProperties(rho, mu);

        // Configure SIMPLE parameters
        simpleSolver.setRelaxationFactors(alphaU, alphaP);
        simpleSolver.setConvergenceTolerance(convergenceTolerance);
        simpleSolver.setMaxIterations(maxIterations);

        std::cout << "SIMPLE parameters:" << std::endl;
        std::cout << "  Max iterations: " << maxIterations << std::endl;
        std::cout << "  Convergence tolerance: " << convergenceTolerance << std::endl;
        std::cout << "  Velocity relaxation: " << alphaU << std::endl;
        std::cout << "  Pressure relaxation: " << alphaP << std::endl;

        // Configure field constraints
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
        const ScalarField* mu_t_field = simpleSolver.getTurbulentViscosity();
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
        if (turbulenceEnabled && k_field && omega_field && mu_t_field && wallDist_field) {
            scalarFieldsToVtk["k"] = k_field;
            scalarFieldsToVtk["omega"] = omega_field;
            scalarFieldsToVtk["mu_t"] = mu_t_field;
            scalarFieldsToVtk["wallDistance"] = wallDist_field;
        }

        // Prepare vector fields for export
        std::map<std::string, const VectorField*> vectorFieldsToVtk;
        vectorFieldsToVtk["velocity"] = &velocity;

        // Write VTK file
        VtkWriter::writeVtkFile(
            vtkOutputFilename,
            allNodes,
            allFaces,
            scalarFieldsToVtk
        );

        std::cout << "Results exported to: " << vtkOutputFilename << std::endl;

        /**********************************************************************
         * -------------------------- 8. FINALIZE ------------------------------
         *********************************************************************/

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

        std::cout << "\n--- Simulation Complete ---" << std::endl;
        std::cout << "Total execution time: " << duration.count() << " seconds" << std::endl;
        std::cout << "Configuration file used: " << configFile << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}