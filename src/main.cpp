#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <map>
#include <algorithm>

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
#include "MatrixConstructor.h"
#include "LinearSolvers.h"

// I/O
#include "MeshReader.h"
#include "VtkWriter.h"

int main() {
    std::cout << "--- Welcome to the FVM Solver ---" << std::endl;
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
        std::cout << "\n--- 2. Setting up Problem ---" << std::endl;

        // --- Physical Properties ---
        const Scalar rho = 1.0;         // Density (kg/m^3)
        const Scalar k = 1.0;           // Thermal conductivity (W/m.K)
        const Scalar cp = 434.0;        // Specific heat (J/kg.K)
        const Scalar Gamma = 0.1;       // Diffusivity coefficient (m^2/s)

        // --- Boundary Conditions ---
        BoundaryConditions bcManager;
        for (const auto& patch : allBoundaryPatches) {
            bcManager.addPatch(patch);
        }

        const std::string phi_field = "phi";
        bcManager.setFixedValue("inlet", phi_field, 1.0);
        bcManager.setFixedValue("outlet", phi_field, 0.0);
        bcManager.printSummary(false);

        // --- Time Control ---
        const Scalar time_total = 1;  // Total simulation time
        const Scalar dt = 0;          // Time step
        int time_step_count = 0;
        const int max_time_steps = static_cast<int>(time_total / dt);

        // --- Discretization Scheme Selection ---
        CentralDifferenceScheme convScheme;  // Use concrete implementation instead of abstract class

        // =========================================================================
        // --- 3. FIELD INITIALIZATION ---
        // =========================================================================
        std::cout << "\n--- 3. Initializing Fields ---" << std::endl;
        ScalarField phi(phi_field, allCells.size(), 0.0); 
        ScalarField phi_old = phi; // Store previous time step values
        VelocityField U("U", allCells.size(), Vector(0, 0, 0.1));

        // =========================================================================
        // --- 4. SOLVER SETUP ---
        // =========================================================================
        TimeScheme timeScheme = TimeScheme::Steady;
        GradientScheme gradScheme;
        MatrixConstructor matrixBuilder(allFaces, allCells, bcManager, gradScheme);

        // =========================================================================
        // --- 5. SOLVER LOOP ---
        // =========================================================================
        if (timeScheme == TimeScheme::Steady) {
            std::cout << "\n--- 5. Starting Steady Solver Loop ---" << std::endl;
            matrixBuilder.constructScalarTransportMatrix(phi_field, phi, phi_old, U, rho, Gamma, TimeScheme::Steady, 0.0, 1.0, convScheme);

            const auto& A_phi = matrixBuilder.getMatrixA();
            const auto& b_phi = matrixBuilder.getVectorB();

            std::cout << "    " << phi_field << " Matrix A (" << A_phi.rows() << "x" << A_phi.cols()
                      << "), Non-zeros: " << A_phi.nonZeros() << std::endl;

            if (A_phi.rows() > 0 && A_phi.rows() <= 200) { 
                std::cout << "\n    --- Coefficients of Matrix A (for Phi) ---" << std::endl;
                // Iterate through all non-zero elements of the sparse matrix
                for (int k_outer = 0; k_outer < A_phi.outerSize(); ++k_outer) {
                    for (Eigen::SparseMatrix<Scalar>::InnerIterator it(A_phi, k_outer); it; ++it) {
                        std::cout << "      A_phi(" << it.row() << ", " << it.col() << ") = " << it.value() << std::endl;
                    }
                }
                std::cout << "\n    --- Coefficients of Vector b (for Phi) ---" << std::endl;
                // Iterate through the dense vector
                for (int i = 0; i < b_phi.size(); ++i) {
                    std::cout << "      b_phi(" << i << ") = " << b_phi(i) << std::endl;
                }
                std::cout << "    --- End of A_phi and b_phi coefficients ---" << std::endl;
            }

            if (A_phi.rows() > 0) {
                Eigen::Matrix<Scalar, Eigen::Dynamic, 1> x_phi(A_phi.rows());
                for (size_t i = 0; i < phi.internalField.size() && i < (size_t)x_phi.size(); ++i) {
                    x_phi(i) = phi.internalField[i];
                }
                bool converged = LinearSolvers::BiCGSTAB(x_phi, A_phi, b_phi, 1e-9, 1000, phi_field);
                if (converged) {
                    for (size_t i = 0; i < phi.internalField.size() && i < (size_t)x_phi.size(); ++i) {
                        phi.internalField[i] = x_phi(i);
                    }
                }
            }
        } else {
            std::cout << "\n--- 5. Starting Transient Solver Loop ---" << std::endl;

            std::cout << "Total Time: " << time_total << "s, Time Step (dt): " << dt << "s" << std::endl;
        
            const Scalar theta = 1.0; // Crank-Nicolson scheme

            for (Scalar time = 0; time < time_total; time += dt) {
                time_step_count++;
                std::cout << "\n--- Time Step: " << time_step_count << ", Time: " << time + dt << "s ---" << std::endl;

                // Store current phi in phi_old before solving for the new phi
                phi_old = phi;

                // Construct the matrix system for phi(n+1)
                matrixBuilder.constructScalarTransportMatrix(
                    phi_field,
                    phi,                      // Current phi field (provides initial guess)
                    phi_old,                  // phi field from previous time step
                    U,                      // Velocity (zero for this problem)
                    rho,                    // Density
                    Gamma,                  // Diffusivity
                    TimeScheme::Steady,
                    dt,                     // Time step
                    theta,
                    convScheme
                );

                // Get the matrix (A) and vector (b)
                const auto& A_phi = matrixBuilder.getMatrixA();
                const auto& b_phi = matrixBuilder.getVectorB();

                std::cout << "    " << phi_field << " Matrix A (" << A_phi.rows() << "x" << A_phi.cols()
                        << "), Non-zeros: " << A_phi.nonZeros() << std::endl;

                if (A_phi.rows() > 0 && A_phi.rows() <= 200) { 
                    std::cout << "\n    --- Coefficients of Matrix A (for phi) ---" << std::endl;
                    // Iterate through all non-zero elements of the sparse matrix
                    for (int k_outer = 0; k_outer < A_phi.outerSize(); ++k_outer) {
                        for (Eigen::SparseMatrix<Scalar>::InnerIterator it(A_phi, k_outer); it; ++it) {
                            std::cout << "      A_phi(" << it.row() << ", " << it.col() << ") = " << it.value() << std::endl;
                        }
                    }
                    std::cout << "\n    --- Coefficients of Vector b (for phi) ---" << std::endl;
                    // Iterate through the dense vector
                    for (int i = 0; i < b_phi.size(); ++i) {
                        std::cout << "      b_phi(" << i << ") = " << b_phi(i) << std::endl;
                    }
                    std::cout << "    --- End of A_phi and b_phi coefficients ---" << std::endl;
                }

                if (A_phi.rows() > 0) {
                    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> x_phi(A_phi.rows());

                    for(size_t i = 0; i < phi.size(); ++i) {
                        if (i < (size_t)x_phi.size()) x_phi(i) = phi[i];
                    }

                    // 3. Call the solver with the correct Eigen vector type.
                    bool converged = LinearSolvers::BiCGSTAB(
                        x_phi,                // Pass the Eigen vector `x_phi`
                        A_phi,
                        b_phi,
                        1e-9,
                        1000,
                        phi_field
                    );

                    // 4. If the solve was successful, copy the results from the
                    //    Eigen vector back into your `phi` field.
                    if (converged) {
                        for (size_t i = 0; i < (size_t)x_phi.size(); ++i) {
                            if (i < phi.size()) phi[i] = x_phi(i);
                        }
                    } else {
                        std::cerr << "Warning: Solver failed to converge at this time step." << std::endl;
                    }
                }
            }
        }

        std::cout << "\n--- Solver Loop Finished ---" << std::endl;

        // =========================================================================
        // --- 6. EXPORT RESULTS ---
        // =========================================================================
        std::cout << "\n--- 6. Exporting Results to VTK ---" << std::endl;
        std::string vtkOutputFilename = "../output_files/steady_convectionDiffusion.vtk";

        std::map<std::string, const ScalarField*> scalarFieldsToVtk;
        scalarFieldsToVtk["phi"] = &phi;

        VtkWriter::writeVtkFile(vtkOutputFilename, allNodes, allFaces, allCells, scalarFieldsToVtk);
        std::cout << "Final results written to " << vtkOutputFilename << std::endl;


    } catch (const std::exception& e) {
        std::cerr << "\n*** A critical error occurred in main: " << e.what() << " ***" << std::endl;
        return 1; // Indicate failure
    }

    std::cout << "\n--- Simulation Complete ---" << std::endl;
    return 0;
}