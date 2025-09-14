/******************************************************************************
 * @file main.cpp
 * @brief Main entry point for the 3D incompressible CFD solver
 * 
 * This file contains the main function that orchestrates the entire CFD
 * simulation process using the SIMPLE algorithm. It handles mesh reading, 
 * boundary condition setup, solver configuration, solution computation, 
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
#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "CellData.h"
#include "FaceData.h"

// Boundary conditions
#include "BoundaryPatch.h"
#include "BoundaryData.h"
#include "BoundaryConditions.h"

// Numerics and Solver
#include "GradientScheme.h"
#include "ConvectionScheme.h"
#include "Matrix.h"
#include "LinearSolvers.h"
#include "SIMPLE.h"
#include "Constraint.h"

// I/O
#include "MeshReader.h"
#include "VtkWriter.h"
#include "checkMesh.h"


/**
 * @brief Main function - Entry point for the CFD simulation
 * 
 * This function coordinates the complete CFD simulation workflow:
 * 1. Mesh reading and geometric property calculation
 * 2. Boundary condition setup for velocity and pressure fields  
 * 3. SIMPLE solver initialization and configuration
 * 4. Steady-state flow solution computation
 * 5. Solution field extraction and statistics
 * 6. Post-processing and VTK export for visualization
 * 7. Solution validation and mass conservation checking
 * 8. Performance timing and summary
 * 
 * The solver uses the SIMPLE algorithm to solve the momentum
 * and pressure- correction equations on unstructured 3D meshes. 
 * Optional turbulence modeling is available via the k-omega SST model.
 * 
 * @return 0 on successful execution, 1 on error
 */
int main() 
{
    // Start timing the total execution
    auto start_time = std::chrono::high_resolution_clock::now();    
    
    std::cout   << "--- Welcome to the CFD Solver with SIMPLE Algorithm ---"
                << std::endl;

    std::cout   << "Running with precision: " << SCALAR_MODE << std::endl;

    std::cout   << std::fixed << std::setprecision(6);

    try 
    {
        /**********************************************************************
         * -------------------------- 1. MESH SETUP ---------------------------
         *********************************************************************/

        std::cout << "\n--- 1. Reading and Preparing Mesh ---" << std::endl;

        /// @brief Container for all mesh nodes/vertices
        std::vector<Vector> allNodes;
        /// @brief Container for all mesh faces (internal + boundary)
        std::vector<Face> allFaces;
        /// @brief Container for all mesh cells/elements
        std::vector<Cell> allCells;
        /// @brief Container for boundary patch definitions
        std::vector<BoundaryPatch> allBoundaryPatches;

        std::string meshFilePath = "../inputFiles/pipe_320k.msh";

        readMshFile
        (
            meshFilePath,
            allNodes,
            allFaces,
            allCells,
            allBoundaryPatches
        );

        std::cout   << "Mesh Loaded: " << allNodes.size() << " nodes, "
                    << allFaces.size() << " faces, " << allCells.size()
                    << " cells." << std::endl;

        for (auto& face : allFaces)
        {
            face.calculateGeometricProperties(allNodes);
            // std::cout << face << std::endl;   // For debugging
        }
        std::cout << "Geometric properties calculated for faces." << std::endl;

        for (auto& cell : allCells)
        {
            cell.calculateGeometricProperties(allFaces);
            // std::cout << cell << std::endl;   // For debugging
        }
        std::cout << "Geometric properties calculated for cells." << std::endl;

        for (auto& face : allFaces)
        {
            face.calculateDistanceProperties(allCells);
        }
        std::cout << "Distance properties calculated for faces." << std::endl;

        // Check mesh min/max area and volume 
        checkMesh(allFaces, allCells);

        /**********************************************************************
         * ------------------------ 2. PROBLEM SETUP --------------------------
         *********************************************************************/

        std::cout << "\n--- 2. Setting up CFD Problem ---" << std::endl;

        // --- Physical Properties ---
        const Scalar rho = 1.225;               // Density (kg/m^3)
        const Scalar mu = 1.7894e-5;            // Dynamic viscosity (Pa·s)

        // --- Boundary Conditions ---
        BoundaryConditions bcManager;
        for (const auto& patch : allBoundaryPatches) 
        {
            bcManager.addPatch(patch);
        }

        const std::string U_field = "U";
        const std::string p_field = "p";
        
        // Inlet: Fixed velocity (low Re for stability)
        bcManager.setFixedValue("inlet", U_field, (0.0, 0.0, -0.1));
        bcManager.setZeroGradient("inlet", p_field);
        
        // Outlet: Fixed pressure, zero gradient velocity
        bcManager.setFixedValue("outlet", p_field, 0.0);
        bcManager.setZeroGradient("outlet", U_field);
        
        // Cylinder wall: no-slip for velocity, zero gradient for pressure
        bcManager.setNoSlip("cylinder", U_field);
        bcManager.setZeroGradient("cylinder", p_field);

        // Sphere wall: no-slip for velocity, zero gradient for pressure
        bcManager.setNoSlip("sphere", U_field);
        bcManager.setZeroGradient("sphere", p_field);
        
        // Symmetry boundaries  
        bcManager.setZeroGradient("symmetry1", U_field);
        bcManager.setZeroGradient("symmetry2", U_field);
        bcManager.setZeroGradient("symmetry3", U_field);
        bcManager.setZeroGradient("symmetry4", U_field);
        bcManager.setZeroGradient("symmetry1", p_field);
        bcManager.setZeroGradient("symmetry2", p_field);
        bcManager.setZeroGradient("symmetry3", p_field);
        bcManager.setZeroGradient("symmetry4", p_field);
        
        // Walls: no-slip velocity and zero gradient for pressure
        bcManager.setNoSlip("wall1", U_field);
        bcManager.setNoSlip("wall2", U_field);
        bcManager.setNoSlip("wall3", U_field);
        bcManager.setNoSlip("wall4", U_field);
        bcManager.setNoSlip("wall5", U_field);
        bcManager.setNoSlip("wall6", U_field);
        bcManager.setZeroGradient("wall1", p_field);
        bcManager.setZeroGradient("wall2", p_field);
        bcManager.setZeroGradient("wall3", p_field);
        bcManager.setZeroGradient("wall4", p_field);
        bcManager.setZeroGradient("wall5", p_field);
        bcManager.setZeroGradient("wall6", p_field);

        // Turbulence boundary conditions helper functions
        auto calculateInletTurbulenceValues = [](Scalar U_inlet, Scalar mu, Scalar rho) -> std::pair<Scalar, Scalar> {
            // Turbulence intensity I = 8% (higher for more significant turbulence)
            const Scalar I = 0.08;

            // Turbulent kinetic energy: k = 1.5 * (I * U)²
            const Scalar k_inlet = 1.5 * (I * U_inlet) * (I * U_inlet);

            // Turbulent length scale: L = 0.07 * D_h (hydraulic diameter)
            // For pipe: D_h = 0.05 m (estimate based on mesh)
            const Scalar D_h = 0.05;
            const Scalar L = 0.07 * D_h;

            // Specific dissipation rate: ω = sqrt(k) / (C_μ^0.25 * L)
            const Scalar C_mu = 0.09;
            const Scalar omega_inlet = std::sqrt(k_inlet) / (std::pow(C_mu, 0.25) * L);

            return {k_inlet, omega_inlet};
        };

        // Calculate turbulence inlet values
        const Scalar U_inlet_mag = 0.1;  // Inlet velocity magnitude
        auto [k_inlet, omega_inlet] = calculateInletTurbulenceValues(U_inlet_mag, mu, rho);

        std::cout << "Turbulence inlet conditions:" << std::endl;
        std::cout << "  k_inlet = " << k_inlet << " m²/s²" << std::endl;
        std::cout << "  omega_inlet = " << omega_inlet << " 1/s" << std::endl;

        // Turbulence boundary conditions
        // Inlet: fixed values for k and omega, zero gradient for phi_wall
        bcManager.setFixedValue("inlet", "k", k_inlet);
        bcManager.setFixedValue("inlet", "omega", omega_inlet);
        bcManager.setZeroGradient("inlet", "phi_wall");

        // Outlet: zero gradient for all turbulence fields
        bcManager.setZeroGradient("outlet", "k");
        bcManager.setZeroGradient("outlet", "omega");
        bcManager.setZeroGradient("outlet", "phi_wall");

        // Walls: k=0, omega wall function (handled by turbulence model), phi_wall=0
        bcManager.setFixedValue("wall1", "k", 0.0);
        bcManager.setFixedValue("wall2", "k", 0.0);
        bcManager.setFixedValue("wall3", "k", 0.0);
        bcManager.setFixedValue("wall4", "k", 0.0);

        // For omega at walls, use zero gradient (wall function applied by turbulence model)
        bcManager.setZeroGradient("wall1", "omega");
        bcManager.setZeroGradient("wall2", "omega");
        bcManager.setZeroGradient("wall3", "omega");
        bcManager.setZeroGradient("wall4", "omega");

        // phi_wall = 0 at walls (zero distance)
        bcManager.setFixedValue("wall1", "phi_wall", 0.0);
        bcManager.setFixedValue("wall2", "phi_wall", 0.0);
        bcManager.setFixedValue("wall3", "phi_wall", 0.0);
        bcManager.setFixedValue("wall4", "phi_wall", 0.0);

        bcManager.printSummary();

        // Inital fields 
        Vector initialVelocity(0.0, 0.0, -0.1);
        Scalar initialPressure = 0.0;


        // --- Discretization Scheme Selection ---
        CentralDifferenceScheme CDS;
        SecondOrderUpwindScheme SOUS;
        UpwindScheme UDS;
        GradientScheme gradScheme;

        /**********************************************************************
         * --------------------- 3. SIMPLE SOLVER SETUP -----------------------
         *********************************************************************/

        std::cout << "\n--- 3. Initializing SIMPLE Solver ---" << std::endl;
        
        SIMPLE simpleSolver(allFaces, allCells, bcManager, gradScheme, UDS);

        // Enable turbulence modeling
        simpleSolver.enableTurbulenceModeling(true);

        // Initialize velocity and pressure fields
        simpleSolver.initialize(initialVelocity, initialPressure);

        simpleSolver.setPhysicalProperties(rho, mu);

        // Configure SIMPLE parameters

        // Under-relaxation: Optimized for stable turbulence convergence
        simpleSolver.setRelaxationFactors(0.3, 0.1);   // Standard relaxation factors
        simpleSolver.setConvergenceTolerance(1e-6);    // Tight tolerance for accuracy
        simpleSolver.setMaxIterations(1);             // Sufficient iterations for convergence
        
        // Configure field constraints: max velocity and pressyre range 
        Constraint* constraintSystem = simpleSolver.getConstraintSystem();
        if (constraintSystem)
        {
            constraintSystem->setVelocityConstraints(0.3);
            constraintSystem->setPressureConstraints(-0.05, 0.05);
            constraintSystem->enableConstraints(false, false);
        }

        /**********************************************************************
         * ------------------- 4. SOLVE STEADY-STATE FLOW ---------------------
         *********************************************************************/
        
        std::cout   << "\n--- 4. Solving Steady-State Flow with SIMPLE ---"
                    << std::endl;

        simpleSolver.solve();

        /**********************************************************************
         * ------------------- 5. EXTRACT SOLUTION FIELDS ---------------------
         *********************************************************************/
        
        std::cout << "\n--- 5. Extracting Solution Fields ---" << std::endl;
        
        const VectorField& velocity = simpleSolver.getVelocity();
        const ScalarField& pressure = simpleSolver.getPressure();
        
        // Extract turbulence fields if available
        const ScalarField* k_field = simpleSolver.getTurbulentKineticEnergy();

        const ScalarField* omega_field =
            simpleSolver.getSpecificDissipationRate();

        const ScalarField* mu_t_field = simpleSolver.getTurbulentViscosity();

        const ScalarField* wallDist_field = simpleSolver.getWallDistance();

        // Extract 3D velocity components
        ScalarField U_x = simpleSolver.getVelocityX();
        ScalarField U_y = simpleSolver.getVelocityY();
        ScalarField U_z = simpleSolver.getVelocityZ();
        
        std::cout   << "Solution extracted:"
                    << std::endl;

        /**********************************************************************
         * ----------------------- 6. POST-PROCESSING -------------------------
         *********************************************************************/

        std::cout << "\n--- 6. Post-Processing Results ---" << std::endl;
        
        // Calculate velocity magnitude
        ScalarField velocityMagnitude("velocityMagnitude", velocity.size());

        for (size_t i = 0; i < velocity.size(); ++i)
        {
            velocityMagnitude[i] = velocity[i].magnitude();
        }
        
        // Print some statistics
        Scalar maxVelocity = 0.0, avgVelocity = 0.0;
        Scalar maxPressure = pressure[0], minPressure = pressure[0];
        
        for (size_t i = 0; i < velocity.size(); ++i)
        {
            Scalar vmag = velocityMagnitude[i];
            maxVelocity = std::max(maxVelocity, vmag);
            avgVelocity += vmag;
            
            maxPressure = std::max(maxPressure, pressure[i]);
            minPressure = std::min(minPressure, pressure[i]);
        }
        avgVelocity /= velocity.size();
        
        std::cout   << "Flow Statistics:"
                    << std::endl;

        std::cout   << "  Max velocity magnitude: "
                    << maxVelocity << " m/s" << std::endl;

        std::cout   << "  Average velocity magnitude: "
                    << avgVelocity << " m/s" << std::endl;

        std::cout   << "  Pressure range: [" << minPressure << ", "
                    << maxPressure << "] Pa" << std::endl;
        
        // 3D velocity component statistics
        Scalar maxUx = 0.0, maxUy = 0.0, maxUz = 0.0;
        Scalar avgUx = 0.0, avgUy = 0.0, avgUz = 0.0;
        
        for (size_t i = 0; i < velocity.size(); ++i)
        {
            maxUx = std::max(maxUx, std::abs(U_x[i]));
            maxUy = std::max(maxUy, std::abs(U_y[i]));
            maxUz = std::max(maxUz, std::abs(U_z[i]));
            avgUx += U_x[i];
            avgUy += U_y[i];
            avgUz += U_z[i];
        }
        avgUx /= velocity.size();
        avgUy /= velocity.size();
        avgUz /= velocity.size();
        
        std::cout   << "3D Velocity Component Statistics:" << std::endl;
        
        std::cout   << "  Max U_x: " << maxUx << " m/s, Avg U_x: "
                    << avgUx << " m/s" << std::endl;

        std::cout   << "  Max U_y: " << maxUy << " m/s, Avg U_y: "
                    << avgUy << " m/s" << std::endl;

        std::cout   << "  Max U_z: " << maxUz << " m/s, Avg U_z: "
                    << avgUz << " m/s" << std::endl;
        
        // Turbulence statistics
        if (k_field && omega_field && mu_t_field)
        {
            Scalar maxK = 0.0, avgK = 0.0;
            Scalar maxOmega = 0.0, avgOmega = 0.0;
            Scalar maxMuT = 0.0, avgMuT = 0.0;
            
            for (size_t i = 0; i < k_field->size(); ++i)
            {
                maxK = std::max(maxK, (*k_field)[i]);
                avgK += (*k_field)[i];
                
                maxOmega = std::max(maxOmega, (*omega_field)[i]);
                avgOmega += (*omega_field)[i];
                
                maxMuT = std::max(maxMuT, (*mu_t_field)[i]);
                avgMuT += (*mu_t_field)[i];
            }
            
            avgK /= k_field->size();
            avgOmega /= omega_field->size();
            avgMuT /= mu_t_field->size();
            
            std::cout   << "Turbulence Statistics:" << std::endl;

            std::cout   << "  Max turbulent kinetic energy: "
                        << maxK << " m²/s²" << std::endl;

            std::cout   << "  Average turbulent kinetic energy: "
                        << avgK << " m²/s²" << std::endl;

            std::cout   << "  Max specific dissipation rate: "
                        << maxOmega << " 1/s" << std::endl;

            std::cout   << "  Average specific dissipation rate: "
                        << avgOmega << " 1/s" << std::endl;

            std::cout   << "  Max turbulent viscosity: "
                        << maxMuT << " Pa·s" << std::endl;

            std::cout   << "  Average turbulent viscosity: "
                        << avgMuT << " Pa·s" << std::endl;

            std::cout   << "  Turbulent viscosity ratio (μₜ/μ): "
                        << avgMuT / mu << std::endl;
        }

        /**********************************************************************
         * ------------------------ 7. EXPORT RESULTS -------------------------
         *********************************************************************/

        std::cout << "\n--- 7. Exporting Results to VTK ---" << std::endl;
        // Create output filename for steady-state solution
        std::string vtkOutputFilename = 
            "../outputFiles/pipe_320k.vtp";

        // Prepare scalar fields for export
        std::map<std::string, const ScalarField*> scalarFieldsToVtk;
        scalarFieldsToVtk["pressure"] = &pressure;
        scalarFieldsToVtk["velocityMagnitude"] = &velocityMagnitude;
        scalarFieldsToVtk["U_x"] = &U_x;
        scalarFieldsToVtk["U_y"] = &U_y;
        scalarFieldsToVtk["U_z"] = &U_z;
        
        // Add turbulence fields if available
        if (k_field && omega_field && mu_t_field && wallDist_field) 
        {
            scalarFieldsToVtk["k"] = k_field;
            scalarFieldsToVtk["omega"] = omega_field;
            scalarFieldsToVtk["mu_t"] = mu_t_field;
            scalarFieldsToVtk["wallDistance"] = wallDist_field;
            
            // Calculate additional turbulence quantities for visualization
            ScalarField turbulentIntensity
            (
                "turbulentIntensity", 
                velocity.size()
            );

            ScalarField turbulentViscosityRatio
            (
                "turbulentViscosityRatio",
                velocity.size()
            );

            ScalarField yPlus("yPlus", velocity.size());
            
            for (size_t i = 0; i < velocity.size(); ++i)
            {
                // Turbulent intensity: I = √(2k/3) / |U|
                Scalar U_mag = velocity[i].magnitude();
                if (U_mag > 1e-10)
                {
                    turbulentIntensity[i] =
                        std::sqrt(2.0 * (*k_field)[i] / 3.0) / U_mag;
                }
                else
                {
                    turbulentIntensity[i] = 0.0;
                }
                
                // Turbulent viscosity ratio: μₜ/μ
                turbulentViscosityRatio[i] = (*mu_t_field)[i] / mu;
                
                // y+ approximation: y+ ≈ ρ * √(τ_wall/ρ) * y / μ
                Scalar y = (*wallDist_field)[i];
                Scalar u_tau_approx = std::sqrt((*k_field)[i]);
                yPlus[i] = rho * u_tau_approx * y / mu;
            }
            
            scalarFieldsToVtk["turbulentIntensity"] = &turbulentIntensity;

            scalarFieldsToVtk["turbulentViscosityRatio"] = 
                &turbulentViscosityRatio;

            scalarFieldsToVtk["yPlus"] = &yPlus;
            
            std::cout   << "  Turbulence fields included in VTK export"
                        << std::endl;
        }

        VtkWriter::writeVtkFile
        (
            vtkOutputFilename,
            allNodes,
            allFaces,
            scalarFieldsToVtk
        );

        std::cout   << "Flow solution written to " 
                    << vtkOutputFilename << std::endl;

        /**********************************************************************
         * ------------------- 8. VALIDATION AND DIAGNOSTICS ------------------
         *********************************************************************/
        
        std::cout << "\n--- 8. Solution Validation ---" << std::endl;
        
        // Check mass conservation
        Scalar totalMassImbalance = 0.0;
        for (size_t i = 0; i < allCells.size(); ++i)
        {
            Scalar cellImbalance = 0.0;
            for (size_t j = 0; j < allCells[i].faceIndices().size(); ++j)
            {
                size_t faceIdx = allCells[i].faceIndices()[j];

                int sign = allCells[i].faceSigns()[j];

                cellImbalance += 
                    sign * simpleSolver.getRhieChowFlowRate()[faceIdx];
            }
            totalMassImbalance += std::abs(cellImbalance);
        }
        
        Scalar avgMassImbalance = totalMassImbalance / allCells.size();

        std::cout   << "Mass Conservation Check:" << std::endl;

        std::cout   << "  Average mass imbalance per cell: "
                    << avgMassImbalance << " kg/s" << std::endl;
        
        if (avgMassImbalance < 1e-8)
        {
            std::cout   << "  Mass conservation satisfied" << std::endl;
        }
        else if (avgMassImbalance < 1e-6)
        {
            std::cout   << "  Mass conservation acceptable" << std::endl;
        }
        else
        {
            std::cout   << "  Mass conservation violated - check solution"
                        << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr   << "\n*** A critical error occurred in main: "
                    << e.what() << " ***" << std::endl;
        return 1;
    }

    // Calculate and print total execution time
    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>
        (
            end_time - start_time
        );
    
    std::cout << "\n--- CFD Simulation Complete ---" << std::endl;
    std::cout << "\n=== EXECUTION TIME SUMMARY ===" << std::endl;
    
    // Format time in hours:minutes:seconds
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);

    auto minutes = 
        std::chrono::duration_cast<std::chrono::minutes>(duration - hours);

    auto seconds = 
        std::chrono::duration_cast<std::chrono::seconds>
        (
            duration - hours - minutes
        );
    
    if (hours.count() > 0)
    {
        std::cout   << "Total execution time: " << hours.count() << "h "
                    << minutes.count() << "m " << seconds.count() << "s"
                    << std::endl;
    }
    else if (minutes.count() > 0)
    {
        std::cout   << "Total execution time: " << minutes.count() << "m " 
                    << seconds.count() << "s" << std::endl;
    }
    else
    {
        std::cout   << "Total execution time: " << seconds.count() << "s ("
                    << std::endl;
    }
    
    return 0;
}