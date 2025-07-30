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
#include "MatrixConstructor.h"
#include "LinearSolvers.h"
#include "SIMPLE.h"

// I/O
#include "MeshReader.h"
#include "VtkWriter.h"

int main() {
    // Start timing the total execution
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "--- Welcome to the CFD Solver with SIMPLE Algorithm ---" << std::endl;
    std::cout << "Running with precision: " << SCALAR_MODE << std::endl;

    std::cout << std::fixed << std::setprecision(6);

    try {
        // =========================================================================
        // --- 1. MESH SETUP ---
        // =========================================================================
        std::cout << "\n--- 1. Reading and Preparing Mesh ---" << std::endl;
        std::vector<Vector> allNodes;
        std::vector<Face> allFaces;
        std::vector<Cell> allCells;
        std::vector<BoundaryPatch> allBoundaryPatches;

        std::string meshFilePath = "../input_files/rod_convection.msh";
        readMshFile(meshFilePath, allNodes, allFaces, allCells, allBoundaryPatches);
        std::cout << "Mesh Loaded: " << allNodes.size() << " nodes, "
                  << allFaces.size() << " faces, " << allCells.size() << " cells." << std::endl;

        for (auto& face : allFaces) {
            face.calculateGeometricProperties(allNodes);
            //std::cout << face << std::endl;   // For debugging
        }

        for (auto& cell : allCells) {
            cell.calculateGeometricProperties(allFaces);
            //std::cout << cell << std::endl;   // For debugging
        }
        std::cout << "Geometric properties calculated." << std::endl;

        // =========================================================================
        // --- 2. PROBLEM SETUP ---
        // =========================================================================
        std::cout << "\n--- 2. Setting up CFD Problem ---" << std::endl;

        // --- Physical Properties ---
        const Scalar rho = 1.25;           // Density (kg/m^3)
        const Scalar mu = 1.8e-5;          // Dynamic viscosity (Pa·s)

        // --- Boundary Conditions ---
        BoundaryConditions bcManager;
        for (const auto& patch : allBoundaryPatches) {
            bcManager.addPatch(patch);
        }

        // Set velocity boundary conditions
        const std::string U_field = "U";
        const std::string p_field = "p";
        
        // Inlet: Fixed velocity
        bcManager.setFixedValue("inlet", U_field, Vector(0.1, 0.0, 0.0));  // 0.1 m/s inlet velocity
        bcManager.setFixedGradient("inlet", p_field, 0.0);  // Zero gradient for pressure at inlet
        
        // Outlet: Fixed pressure, zero gradient velocity
        bcManager.setFixedValue("outlet", p_field, 0.0);   // Reference pressure at outlet
        bcManager.setFixedGradient("outlet", U_field, Vector(0.0, 0.0, 0.0));  // Zero gradient for velocity
        
        // Walls: No-slip condition for velocity, zero gradient for pressure
        bcManager.setFixedValue("wall", U_field, Vector(0.0, 0.0, 0.0));  // No-slip wall
        bcManager.setFixedGradient("wall", p_field, 0.0);  // Zero gradient for pressure at wall
        
        bcManager.printSummary(false);

        // --- Discretization Scheme Selection ---
        CentralDifferenceScheme convScheme;  // Use central difference scheme
        GradientScheme gradScheme;                  // Gradient calculation scheme

        // =========================================================================
        // --- 3. SIMPLE SOLVER SETUP ---
        // =========================================================================
        std::cout << "\n--- 3. Initializing SIMPLE Solver with k-omega SST Turbulence ---" << std::endl;
        
        SIMPLE simpleSolver(allFaces, allCells, bcManager, gradScheme, convScheme);
        
        // Configure SIMPLE parameters
        simpleSolver.setRelaxationFactors(0.7, 0.3);  // Under-relaxation: U=0.7, p=0.3
        simpleSolver.setConvergenceTolerance(1e-6);   // Convergence tolerance
        simpleSolver.setMaxIterations(500);           // Maximum iterations
        
        // Enable k-omega SST turbulence modeling
        simpleSolver.enableTurbulenceModeling(true);

        // =========================================================================
        // --- 4. SOLVE PRESSURE-VELOCITY COUPLING WITH TURBULENCE ---
        // =========================================================================
        std::cout << "\n--- 4. Solving Incompressible Turbulent Flow with SIMPLE + k-omega SST ---" << std::endl;
        
        // Solve the coupled pressure-velocity-turbulence system
        simpleSolver.solve();

        // =========================================================================
        // --- 5. EXTRACT SOLUTION FIELDS ---
        // =========================================================================
        std::cout << "\n--- 5. Extracting Solution Fields ---" << std::endl;
        
        const VectorField& velocity = simpleSolver.getVelocity();
        const ScalarField& pressure = simpleSolver.getPressure();
        const FaceFluxField& massFlux = simpleSolver.getMassFlux();
        
        // Extract turbulence fields if available
        const ScalarField* k_field = simpleSolver.getTurbulentKineticEnergy();
        const ScalarField* omega_field = simpleSolver.getSpecificDissipationRate();
        const ScalarField* mu_t_field = simpleSolver.getTurbulentViscosity();
        const ScalarField* wallDist_field = simpleSolver.getWallDistance();
        
        std::cout << "Solution extracted:" << std::endl;
        std::cout << "  Velocity field size: " << velocity.size() << std::endl;
        std::cout << "  Pressure field size: " << pressure.size() << std::endl;
        std::cout << "  Mass flux field size: " << massFlux.size() << std::endl;
        
        if (k_field && omega_field && mu_t_field) {
            std::cout << "  Turbulence fields available:" << std::endl;
            std::cout << "    k field size: " << k_field->size() << std::endl;
            std::cout << "    omega field size: " << omega_field->size() << std::endl;
            std::cout << "    mu_t field size: " << mu_t_field->size() << std::endl;
            std::cout << "    Wall distance field size: " << wallDist_field->size() << std::endl;
        }

        // =========================================================================
        // --- 6. POST-PROCESSING ---
        // =========================================================================
        std::cout << "\n--- 6. Post-Processing Results ---" << std::endl;
        
        // Calculate velocity magnitude
        ScalarField velocityMagnitude("velocityMagnitude", velocity.size());
        for (size_t i = 0; i < velocity.size(); ++i) {
            velocityMagnitude[i] = velocity[i].magnitude();
        }
        
        // Calculate pressure coefficient (Cp)
        ScalarField pressureCoeff("pressureCoeff", pressure.size());
        Scalar pRef = 0.0;  // Reference pressure
        Scalar uRef = 0.1;  // Reference velocity (inlet velocity)
        Scalar dynamicPressure = 0.5 * rho * uRef * uRef;
        
        for (size_t i = 0; i < pressure.size(); ++i) {
            pressureCoeff[i] = (pressure[i] - pRef) / dynamicPressure;
        }
        
        // Print some statistics
        Scalar maxVelocity = 0.0, avgVelocity = 0.0;
        Scalar maxPressure = pressure[0], minPressure = pressure[0];
        
        for (size_t i = 0; i < velocity.size(); ++i) {
            Scalar vmag = velocityMagnitude[i];
            maxVelocity = std::max(maxVelocity, vmag);
            avgVelocity += vmag;
            
            maxPressure = std::max(maxPressure, pressure[i]);
            minPressure = std::min(minPressure, pressure[i]);
        }
        avgVelocity /= velocity.size();
        
        std::cout << "Flow Statistics:" << std::endl;
        std::cout << "  Max velocity magnitude: " << maxVelocity << " m/s" << std::endl;
        std::cout << "  Average velocity magnitude: " << avgVelocity << " m/s" << std::endl;
        std::cout << "  Pressure range: [" << minPressure << ", " << maxPressure << "] Pa" << std::endl;
        
        // Turbulence statistics
        if (k_field && omega_field && mu_t_field) {
            Scalar maxK = 0.0, avgK = 0.0;
            Scalar maxOmega = 0.0, avgOmega = 0.0;
            Scalar maxMuT = 0.0, avgMuT = 0.0;
            
            for (size_t i = 0; i < k_field->size(); ++i) {
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
            
            std::cout << "Turbulence Statistics:" << std::endl;
            std::cout << "  Max turbulent kinetic energy: " << maxK << " m²/s²" << std::endl;
            std::cout << "  Average turbulent kinetic energy: " << avgK << " m²/s²" << std::endl;
            std::cout << "  Max specific dissipation rate: " << maxOmega << " 1/s" << std::endl;
            std::cout << "  Average specific dissipation rate: " << avgOmega << " 1/s" << std::endl;
            std::cout << "  Max turbulent viscosity: " << maxMuT << " Pa·s" << std::endl;
            std::cout << "  Average turbulent viscosity: " << avgMuT << " Pa·s" << std::endl;
            std::cout << "  Turbulent viscosity ratio (μₜ/μ): " << avgMuT / mu << std::endl;
        }

        // =========================================================================
        // --- 7. EXPORT RESULTS ---
        // =========================================================================
        std::cout << "\n--- 7. Exporting Results to VTK ---" << std::endl;
        std::string vtkOutputFilename = "../output_files/komega_sst_flow_solution.vtk";

        // Prepare scalar fields for export
        std::map<std::string, const ScalarField*> scalarFieldsToVtk;
        scalarFieldsToVtk["pressure"] = &pressure;
        scalarFieldsToVtk["velocityMagnitude"] = &velocityMagnitude;
        scalarFieldsToVtk["pressureCoeff"] = &pressureCoeff;

        // Prepare vector fields for export (convert velocity to separate components)
        ScalarField U_x("U_x", velocity.size());
        ScalarField U_y("U_y", velocity.size());
        ScalarField U_z("U_z", velocity.size());
        
        for (size_t i = 0; i < velocity.size(); ++i) {
            U_x[i] = velocity[i].x;
            U_y[i] = velocity[i].y;
            U_z[i] = velocity[i].z;
        }
        
        scalarFieldsToVtk["U_x"] = &U_x;
        scalarFieldsToVtk["U_y"] = &U_y;
        scalarFieldsToVtk["U_z"] = &U_z;
        
        // Add turbulence fields to export if available
        if (k_field && omega_field && mu_t_field && wallDist_field) {
            scalarFieldsToVtk["k"] = k_field;
            scalarFieldsToVtk["omega"] = omega_field;
            scalarFieldsToVtk["mu_t"] = mu_t_field;
            scalarFieldsToVtk["wallDistance"] = wallDist_field;
            
            // Calculate additional turbulence quantities for visualization
            ScalarField turbulentIntensity("turbulentIntensity", velocity.size());
            ScalarField turbulentViscosityRatio("turbulentViscosityRatio", velocity.size());
            ScalarField yPlus("yPlus", velocity.size());
            
            for (size_t i = 0; i < velocity.size(); ++i) {
                // Turbulent intensity: I = √(2k/3) / |U|
                Scalar U_mag = velocity[i].magnitude();
                if (U_mag > 1e-10) {
                    turbulentIntensity[i] = std::sqrt(2.0 * (*k_field)[i] / 3.0) / U_mag;
                } else {
                    turbulentIntensity[i] = 0.0;
                }
                
                // Turbulent viscosity ratio: μₜ/μ
                turbulentViscosityRatio[i] = (*mu_t_field)[i] / mu;
                
                // y+ approximation: y+ ≈ ρ * √(τ_wall/ρ) * y / μ
                Scalar y = (*wallDist_field)[i];
                Scalar u_tau_approx = std::sqrt((*k_field)[i]);  // Rough approximation
                yPlus[i] = rho * u_tau_approx * y / mu;
            }
            
            scalarFieldsToVtk["turbulentIntensity"] = &turbulentIntensity;
            scalarFieldsToVtk["turbulentViscosityRatio"] = &turbulentViscosityRatio;
            scalarFieldsToVtk["yPlus"] = &yPlus;
            
            std::cout << "  Turbulence fields included in VTK export" << std::endl;
        }

        VtkWriter::writeVtkFile(vtkOutputFilename, allNodes, allFaces, allCells, scalarFieldsToVtk);
        std::cout << "k-omega SST turbulent flow solution written to " << vtkOutputFilename << std::endl;

        // =========================================================================
        // --- 8. VALIDATION AND DIAGNOSTICS ---
        // =========================================================================
        std::cout << "\n--- 8. Solution Validation ---" << std::endl;
        
        // Check mass conservation
        Scalar totalMassImbalance = 0.0;
        for (size_t i = 0; i < allCells.size(); ++i) {
            Scalar cellImbalance = 0.0;
            for (size_t j = 0; j < allCells[i].faceIndices.size(); ++j) {
                size_t faceIdx = allCells[i].faceIndices[j];
                int sign = allCells[i].faceSigns[j];
                cellImbalance += sign * massFlux[faceIdx];
            }
            totalMassImbalance += std::abs(cellImbalance);
        }
        
        Scalar avgMassImbalance = totalMassImbalance / allCells.size();
        std::cout << "Mass Conservation Check:" << std::endl;
        std::cout << "  Average mass imbalance per cell: " << avgMassImbalance << " kg/s" << std::endl;
        
        if (avgMassImbalance < 1e-8) {
            std::cout << "  ✓ Mass conservation satisfied" << std::endl;
        } else if (avgMassImbalance < 1e-6) {
            std::cout << "  ⚠ Mass conservation acceptable" << std::endl;
        } else {
            std::cout << "  ✗ Mass conservation violated - check solution" << std::endl;
        }

        // Reynolds number calculation
        Scalar characteristicLength = 1.0;  // Assume 1m characteristic length
        Scalar Re = rho * avgVelocity * characteristicLength / mu;
        std::cout << "Flow Parameters:" << std::endl;
        std::cout << "  Reynolds number: " << Re << std::endl;
        
        if (Re < 1) {
            std::cout << "  Flow regime: Stokes flow (very low Re)" << std::endl;
        } else if (Re < 100) {
            std::cout << "  Flow regime: Laminar flow (low Re)" << std::endl;
        } else if (Re < 2300) {
            std::cout << "  Flow regime: Laminar flow" << std::endl;
        } else {
            std::cout << "  Flow regime: Transitional/Turbulent (high Re)" << std::endl;
        }
        
        // Turbulent flow parameters
        if (k_field && omega_field && mu_t_field) {
            Scalar avgK = 0.0, avgOmega = 0.0;
            for (size_t i = 0; i < k_field->size(); ++i) {
                avgK += (*k_field)[i];
                avgOmega += (*omega_field)[i];
            }
            avgK /= k_field->size();
            avgOmega /= omega_field->size();
            
            // Turbulent Reynolds number: Re_t = ρ * k / (μ * ω)
            Scalar Re_t = rho * avgK / (mu * avgOmega + 1e-20);
            
            // Turbulent time scale: T_t = 1/ω
            Scalar T_t = 1.0 / (avgOmega + 1e-20);
            
            // Turbulent length scale: L_t = √k / ω
            Scalar L_t = std::sqrt(avgK) / (avgOmega + 1e-20);
            
            std::cout << "Turbulent Flow Parameters:" << std::endl;
            std::cout << "  Turbulent Reynolds number: " << Re_t << std::endl;
            std::cout << "  Turbulent time scale: " << T_t << " s" << std::endl;
            std::cout << "  Turbulent length scale: " << L_t << " m" << std::endl;
            std::cout << "  Average turbulent intensity: " << std::sqrt(2.0*avgK/3.0)/avgVelocity*100.0 << " %" << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "\n*** A critical error occurred in main: " << e.what() << " ***" << std::endl;
        return 1; // Indicate failure
    }

    // Calculate and print total execution time
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "\n--- CFD Simulation with SIMPLE + k-omega SST Turbulence Complete ---" << std::endl;
    std::cout << "\n=== EXECUTION TIME SUMMARY ===" << std::endl;
    std::cout << "Total execution time: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "Total execution time: " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    // Format time in hours:minutes:seconds for longer runs
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration - hours);
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration - hours - minutes);
    
    if (hours.count() > 0) {
        std::cout << "Total execution time: " << hours.count() << "h " 
                  << minutes.count() << "m " << seconds.count() << "s" << std::endl;
    } else if (minutes.count() > 0) {
        std::cout << "Total execution time: " << minutes.count() << "m " 
                  << seconds.count() << "s" << std::endl;
    }
    
    return 0;
}