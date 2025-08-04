#include "SIMPLE.h"
#include "Scalar.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <map>

SIMPLE::SIMPLE(const std::vector<Face>& faces,
               const std::vector<Cell>& cells, 
               const BoundaryConditions& bc,
               const GradientScheme& gradScheme,
               const ConvectionDiscretization& convScheme)
    : allFaces(faces),
      allCells(cells),
      bcManager(bc),
      gradientScheme(gradScheme),
      convectionScheme(convScheme),
      rho(1.225),
      mu(1.7894e-5),
      alpha_U(0.7),
      alpha_p(0.3),
      maxIterations(500),
      tolerance(1e-6),
      enableTurbulence(false),
      U("U", cells.size(), Vector(0.0, 0.0, 0.0)),
      p("p", cells.size(), 0.0),
      p_prime("p_prime", cells.size(), 0.0),
      U_face("U_face", faces.size(), Vector(0.0, 0.0, 0.0)),
      massFlux("massFlux", faces.size(), 0.0),
      volumeFlux("volumeFlux", faces.size(), 0.0),
      a_U("a_U", cells.size(), 0.0),
      H_U("H_U", cells.size(), Vector(0.0, 0.0, 0.0)),
      grad_P("grad_P", cells.size(), Vector(0.0, 0.0, 0.0))
{
    // Initialize with zero velocity to avoid conflicts with boundary conditions
    initialize(Vector(0.0, 0.0, 0.0), 0.0);
}

void SIMPLE::initialize(const Vector& initialVelocity, Scalar initialPressure) {
    // Initialize matrix constructor
    matrixConstruct = std::make_unique<Matrix>(
        allFaces, allCells, bcManager, gradientScheme);
    
    // Initialize velocity and pressure fields
    for (size_t i = 0; i < allCells.size(); ++i) {
        U[i] = initialVelocity;
        p[i] = initialPressure;                     
    }
    
    // Initialize turbulence model if enabled
    if (enableTurbulence) {
        turbulenceModel = std::make_unique<KOmegaSST>(
            allFaces, allCells, bcManager, gradientScheme);
        turbulenceModel->initialize(U, rho, mu);
        std::cout << "k-omega SST turbulence model initialized." << std::endl;
    }
    
    // Apply initial boundary conditions
    applyVelocityBoundaryConditions();
    applyPressureBoundaryConditions();
    
    std::cout << "SIMPLE algorithm initialized with " << allCells.size() 
              << " cells and " << allFaces.size() << " faces." << std::endl;
}

void SIMPLE::solve() {

    std::cout << "\n=== Starting SIMPLE Loop" << " ===" << std::endl;
    
    int iteration = 0;
    bool converged = false;
    
    // Initial cache build
    matrixConstruct->refreshIterationCaches(p, U, rho, turbulenceModel.get());
    grad_P = matrixConstruct->gradP;
    
    while (!converged && iteration < maxIterations) {
        std::cout << "\n--- SIMPLE Iteration " << iteration + 1 << " ---" << std::endl;

        // Refresh caches (gradients, mdot) for this iteration
        matrixConstruct->refreshIterationCaches(p, U, rho, turbulenceModel.get());
        grad_P = matrixConstruct->gradP;
        
        // Step 1: Solve momentum equations with effective viscosity
        solveMomentumEquations();
            
        // Step 2: Calculate face fluxes using Rhie-Chow interpolation
        calculateRhieChowFaceVelocities();
        calculateMassFluxes();
        
        // Step 3: Solve pressure correction equation
        solvePressureCorrection();
        
        // Step 4: Correct velocities and pressure
        correctVelocity();
        correctPressure();
        
        // Step 5: Solve turbulence equations (k-omega SST sequence)
        if (enableTurbulence && turbulenceModel) {
            std::cout << "  Solving turbulence equations..." << std::endl;
            
            // Solve complete k-omega SST model using cached velocity gradients
            turbulenceModel->solve(U,
                                 matrixConstruct->gradUx,
                                 matrixConstruct->gradUy,
                                 matrixConstruct->gradUz,
                                 rho, mu);
        }
        
        // Step 6: Update boundary conditions
        applyVelocityBoundaryConditions();
        applyPressureBoundaryConditions();
        
        // Step 7: Check convergence
        converged = checkConvergence();
        
        iteration++;
        
        if (iteration % 10 == 0 || converged) {
            std::cout << "Iteration: " << iteration;
            if (converged) std::cout << " - CONVERGED";
            std::cout << std::endl;
        }
    }
    
    if (!converged) {
        std::cout << "WARNING: SIMPLE algorithm did not converge after " 
                  << maxIterations << " iterations." << std::endl;
    } else {
        std::cout << "SIMPLE algorithm converged in " << iteration << " iterations." << std::endl;
    }
}

void SIMPLE::solveMomentumEquations() {
    /**
     * Solve momentum equations with effective viscosity (laminar + turbulent):
     * ∂(ρU)/∂t + ∇·(ρUU) = -∇p + ∇·[(μ + μₜ)∇U] + S
     * 
     * Where μₜ is the turbulent viscosity from k-omega SST model
     */
    
    // Obtain pressure gradient from cache
    grad_P = matrixConstruct->gradP;
    
    // Reset momentum equation diagonal accumulator
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_U[i] = 0.0;
    }
    
    // Get effective viscosity (laminar + turbulent)
    Scalar mu_eff = mu;  // Default to laminar viscosity
    ScalarField mu_effective("mu_eff", allCells.size(), mu);
    
    if (enableTurbulence && turbulenceModel) {
        mu_effective = turbulenceModel->getEffectiveViscosity(mu);
        std::cout << "  Using effective viscosity (μ_lam + μ_t)" << std::endl;
    }
    
    // Extract U.x component and old values for transient
    ScalarField U_x("U_x", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_x[i] = U[i].x;
    }
    
    // Solve U-momentum equation with effective viscosity
    matrixConstruct->constructScalarTransportMatrix(
        "U_x", 
        U_x,  // Current U.x component
        ScalarField("U_x_old", allCells.size()), // Previous time step
        U, // Transport velocity - use actual velocity field
        ScalarField("U_x_source", allCells.size()), // Source term (zero for momentum)
        rho,    // Density
        mu_eff, // Effective viscosity (laminar + turbulent)
        TimeScheme::Steady,
        0.0,    // dt
        1.0, // theta
        convectionScheme
    );
    
    // Old values already extracted above
    
    // Add pressure gradient source term and under-relaxation
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstruct->getVectorB());
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstruct->getMatrixA());
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar cellVolume = allCells[i].volume;
        
        // Add pressure gradient source term: -∇p·V
        b_vector(i) -= grad_P[i].x * cellVolume;
        
        
        // Apply under-relaxation: a_P = a_P/α + (1-α)/α * a_P_old
        Scalar a_P = A_matrix.coeff(i, i);
        Scalar a_P_relaxed = a_P / alpha_U;
        a_U[i] += a_P_relaxed;   // Accumulate for Rhie-Chow (will average later)
        A_matrix.coeffRef(i, i) = a_P_relaxed;
        b_vector(i) += (1.0 - alpha_U) / alpha_U * a_P * U[i].x;
        
    }
    
    // Solve for U.x
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> U_x_solution(allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_x_solution(i) = U[i].x;
    }
    
    bool solved = LinearSolvers::BiCGSTAB(
        U_x_solution, A_matrix, b_vector, 1e-8, 1000, "U_x");
    
    if (solved) {
        for (size_t i = 0; i < allCells.size(); ++i) {
            U[i].x = U_x_solution(i);
        }
    }
    
    // Solve V-momentum equation (y-component)
    ScalarField U_y("U_y", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_y[i] = U[i].y;
    }
    
    matrixConstruct->constructScalarTransportMatrix(
        "U_y", U_y, ScalarField("U_y_old", allCells.size()),
        U, // Transport velocity - use actual velocity field
        ScalarField("U_y_source", allCells.size()), // Source term (zero for momentum)
        rho, mu_eff, 
        TimeScheme::Steady,
        0.0, 
        1.0, 
        convectionScheme);
    
    auto& b_vector_y = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstruct->getVectorB());
    auto& A_matrix_y = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstruct->getMatrixA());
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar cellVolume = allCells[i].volume;
        b_vector_y(i) -= grad_P[i].y * cellVolume;
        
        Scalar a_P = A_matrix_y.coeff(i, i);
        Scalar a_P_relaxed = a_P / alpha_U;
        a_U[i] += a_P_relaxed;
        A_matrix_y.coeffRef(i, i) = a_P_relaxed;
        b_vector_y(i) += (1.0 - alpha_U) / alpha_U * a_P * U[i].y;
    }
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> U_y_solution(allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_y_solution(i) = U[i].y;
    }
    
    solved = LinearSolvers::BiCGSTAB(
        U_y_solution, A_matrix_y, b_vector_y, 1e-8, 1000, "U_y");
    
    if (solved) {
        for (size_t i = 0; i < allCells.size(); ++i) {
            U[i].y = U_y_solution(i);
        }
    }
    
    // Solve W-momentum equation (z-component)
    ScalarField U_z("U_z", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_z[i] = U[i].z;
    }
    
    matrixConstruct->constructScalarTransportMatrix(
        "U_z", U_z, ScalarField("U_z_old", allCells.size()),
        U, // Transport velocity - use actual velocity field
        ScalarField("U_z_source", allCells.size()), // Source term (zero for momentum)
        rho, mu_eff, 
        TimeScheme::Steady,
        0.0, 
        1.0, 
        convectionScheme);
    
    auto& b_vector_z = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstruct->getVectorB());
    auto& A_matrix_z = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstruct->getMatrixA());
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar cellVolume = allCells[i].volume;
        b_vector_z(i) -= grad_P[i].z * cellVolume;
        
        Scalar a_P = A_matrix_z.coeff(i, i);
        Scalar a_P_relaxed = a_P / alpha_U;
        a_U[i] += a_P_relaxed;
        A_matrix_z.coeffRef(i, i) = a_P_relaxed;
        b_vector_z(i) += (1.0 - alpha_U) / alpha_U * a_P * U[i].z;  
    }
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> U_z_solution(allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_z_solution(i) = U[i].z;
    }
    
    solved = LinearSolvers::BiCGSTAB(
        U_z_solution, A_matrix_z, b_vector_z, 1e-8, 1000, "U_z");
    
    if (solved) {
        for (size_t i = 0; i < allCells.size(); ++i) {
            U[i].z = U_z_solution(i);
        }
    }

    // Average the accumulated diagonal from the three momentum components
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_U[i] /= 3.0;
    }
}

void SIMPLE::calculateRhieChowFaceVelocities() {
    // Rhie-Chow interpolation for face velocities
    // U_f = U_f_interpolated + D_f * (∇p_cell_interpolated - ∇p_face_from_neighbors)
    
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx) {
        const Face& face = allFaces[faceIdx];
        
        if (face.isBoundary()) {
            // --- Boundary face handling ---
            size_t P = face.ownerCell;
            Vector U_b(0.0, 0.0, 0.0);

            // Build (or reuse) face-to-patch map
            static std::map<size_t, const BoundaryPatch*> faceToPatchMap;
            if (faceToPatchMap.empty()) {
                for (const auto& patch : bcManager.patches) {
                    for (size_t f = patch.firstFaceIndex; f <= patch.lastFaceIndex; ++f) {
                        faceToPatchMap[f] = &patch;
                    }
                }
            }

            const BoundaryPatch* patch = faceToPatchMap.at(face.id);
            const BoundaryData* bc = bcManager.getFieldBC(patch->patchName, "U");

            if (bc) {
                switch (bc->type) {
                    case BCType::FIXED_VALUE:
                        if (bc->valueType == BCValueType::VECTOR) {
                            U_b = bc->vectorValue;
                        } else {
                            // Treat scalar value as magnitude normal to the face
                            U_b = face.normal * bc->scalarValue;
                        }
                        break;
                    case BCType::NO_SLIP:
                        U_b = Vector(0.0, 0.0, 0.0);
                        break;
                    case BCType::ZERO_GRADIENT:
                        U_b = U[P];
                        break;
                    case BCType::FIXED_GRADIENT: {
                        Vector d_Pf = face.centroid - allCells[P].centroid;
                        if (bc->gradientType == BCValueType::VECTOR) {
                            U_b = U[P] + bc->vectorGradient * d_Pf.magnitude();
                        } else {
                            U_b = U[P] + face.normal * bc->scalarGradient * d_Pf.magnitude();
                        }
                        break;
                    }
                    default:
                        U_b = U[P];
                        break;
                }
            } else {
                // Fallback to owner cell velocity
                U_b = U[P];
            }

            U_face[faceIdx] = U_b;
            continue;
        }
        
        size_t P = face.ownerCell;
        size_t N = face.neighbourCell.value();
        
        // Distance weighting for interpolation
        Vector d_PN = allCells[N].centroid - allCells[P].centroid;
        Vector d_Pf = face.centroid - allCells[P].centroid;
        Vector d_Nf = face.centroid - allCells[N].centroid;
        
        Scalar d_P = d_Pf.magnitude();
        Scalar d_N = d_Nf.magnitude();
        Scalar w_f = d_N / (d_P + d_N);  // Interpolation weight
        
        // Standard interpolated face velocity
        Vector U_f_interpolated = w_f * U[P] + (1.0 - w_f) * U[N];
        
        // Interpolated pressure gradient
        Vector grad_P_f_interpolated = w_f * grad_P[P] + (1.0 - w_f) * grad_P[N];
        
        // Face-based pressure gradient (using neighbor cell pressures)
        Vector gradP_f_face = (p[N] - p[P]) / d_PN.magnitude() * 
                              (d_PN / d_PN.magnitude());
        
        // Interpolated diffusion coefficient D_f
        Scalar a_P_interp = w_f * a_U[P] + (1.0 - w_f) * a_U[N];
        Scalar D_f = face.area / (a_P_interp + 1e-20);  // Avoid division by zero
        
        // Rhie-Chow correction
        Vector correction = D_f * (grad_P_f_interpolated - gradP_f_face);
        
        // Final face velocity
        U_face[faceIdx] = U_f_interpolated - correction;
    }
}

void SIMPLE::calculateMassFluxes() {
    // Calculate mass flux through each face using Rhie-Chow face velocities
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx) {
        const Face& face = allFaces[faceIdx];
        
        // Volume flux: U_f · S_f
        volumeFlux[faceIdx] = dot(U_face[faceIdx], face.normal * face.area);
        
        // Mass flux: ρ * (U_f · S_f)
        massFlux[faceIdx] = rho * volumeFlux[faceIdx];
    }
}

void SIMPLE::solvePressureCorrection() {
    // Build pressure correction equation
    // ∇·(D_f ∇p') = ∇·ρU
    
    matrixConstruct->constructScalarTransportMatrix(
        "p_prime", p_prime, ScalarField("p_prime_old", allCells.size()),
        VectorField("zero_velocity", allCells.size(), Vector(0.0, 0.0, 0.0)),
        ScalarField("p_prime_source", allCells.size()), // Source term (zero for pressure correction)
        0.0,  // No convection for pressure correction
        1.0,  // Diffusion coefficient (will be overridden)
        TimeScheme::Steady, 0.0, 1.0, convectionScheme); // Pressure correction always steady
    
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstruct->getMatrixA());
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstruct->getVectorB());
    
    // Clear the matrix and rebuild for pressure correction
    A_matrix.setZero();
    b_vector.setZero();
    
    std::vector<Eigen::Triplet<Scalar>> triplets;
    
    // Loop over all faces to build pressure correction equation
    for (const auto& face : allFaces) {
        if (!face.geometricPropertiesCalculated) continue;
        
        size_t P = face.ownerCell;
        
        if (face.isBoundary()) {
            // For boundary faces, we need to handle pressure outlets
            // Build face-to-patch map if not already done
            static std::map<size_t, const BoundaryPatch*> faceToPatchMap;
            if (faceToPatchMap.empty()) {
                for (const auto& patch : bcManager.patches) {
                    for (size_t i = patch.firstFaceIndex; i <= patch.lastFaceIndex; ++i) {
                        faceToPatchMap[i] = &patch;
                    }
                }
            }
            
            const BoundaryPatch* patch = faceToPatchMap.at(face.id);
            const BoundaryData* bc = bcManager.getFieldBC(patch->patchName, "p");
            
            if (bc && bc->type == BCType::FIXED_VALUE) {
                // Fixed pressure boundary: p' = 0
                // Add large coefficient to diagonal to enforce p'[P] = 0
                triplets.emplace_back(P, P, 1e10);
                // No contribution to RHS as p' should be zero
            } else {
                // For other boundaries (walls, inlets), treat as zero flux
                // This is implicitly handled by not adding face contribution
            }
        } else {
            // Internal face
            size_t N = face.neighbourCell.value();
            
            // Distance between cell centers
            Vector d_PN = allCells[N].centroid - allCells[P].centroid;
            Scalar d_mag = d_PN.magnitude();
            
            // Pressure correction diffusion coefficient
            Scalar a_P_interp = 0.5 * (a_U[P] + a_U[N]);
            Scalar D_f = face.area / (a_P_interp + 1e-20);
            
            // Matrix coefficients
            triplets.emplace_back(P, P, D_f / d_mag);
            triplets.emplace_back(P, N, -D_f / d_mag);
            triplets.emplace_back(N, N, D_f / d_mag);
            triplets.emplace_back(N, P, -D_f / d_mag);
            
            // Source term: mass imbalance
            b_vector(P) -= massFlux[face.id];
            b_vector(N) += massFlux[face.id];
        }
    }
    
    A_matrix.setFromTriplets(triplets.begin(), triplets.end());
    
    // Fix one pressure correction value to avoid singular matrix
    A_matrix.coeffRef(0, 0) += 1e12;
    b_vector(0) = 0.0;
    
    // Solve pressure correction equation
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> p_prime_solution(allCells.size());
    p_prime_solution.setZero();
    
    bool solved = LinearSolvers::BiCGSTAB(
        p_prime_solution, A_matrix, b_vector, 1e-8, 1000, "p_prime");
    
    if (solved) {
        for (size_t i = 0; i < allCells.size(); ++i) {
            p_prime[i] = p_prime_solution(i);
        }
    }
}

void SIMPLE::correctVelocity() {
    // Velocity correction: U = U* - (V/a_U) * ∇p'
    VectorField gradP_prime = gradientScheme.LeastSquares(p_prime, allCells);
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar cellVolume = allCells[i].volume;
        Vector correction = (cellVolume / (a_U[i] + 1e-20)) * gradP_prime[i];
        U[i] = U[i] - correction;
    }
}

void SIMPLE::correctPressure() {
    // Pressure correction: p = p + α_p * p'
    for (size_t i = 0; i < allCells.size(); ++i) {
        p[i] += alpha_p * p_prime[i];
    }
    
    // Reset pressure correction
    for (size_t i = 0; i < allCells.size(); ++i) {
        p_prime[i] = 0.0;
    }
}

bool SIMPLE::checkConvergence() {
    // Check mass conservation and velocity/pressure residuals
    Scalar massImbalance = calculateMassImbalance();
    Scalar velocityResidual = calculateVelocityResidual();
    Scalar pressureResidual = calculatePressureResidual();
    
    bool converged = (massImbalance < tolerance) && 
                     (velocityResidual < tolerance) && 
                     (pressureResidual < tolerance);
    
    std::cout << "  Residuals - Mass: " << std::scientific << massImbalance 
              << ", Velocity: " << velocityResidual
              << ", Pressure: " << pressureResidual 
              << std::fixed << std::endl;
    
    // Additional convergence monitoring
    if (massImbalance > 1e3 || velocityResidual > 1e3 || pressureResidual > 1e3) {
        std::cout << "  WARNING: Residuals are very large - solution may be diverging!" << std::endl;
        std::cout << "  Consider: " << std::endl;
        std::cout << "    - Reducing relaxation factors (current: U=" << alpha_U << ", p=" << alpha_p << ")" << std::endl;
        std::cout << "    - Checking boundary conditions" << std::endl;
        std::cout << "    - Improving initial guess" << std::endl;
    }
    
    return converged;
}

void SIMPLE::applyVelocityBoundaryConditions() {
    // Build face-to-patch map for boundary lookup
    std::map<size_t, const BoundaryPatch*> faceToPatchMap;
    for (const auto& patch : bcManager.patches) {
        for (size_t i = patch.firstFaceIndex; i <= patch.lastFaceIndex; ++i) {
            faceToPatchMap[i] = &patch;
        }
    }
    
    // Apply velocity boundary conditions
    for (const auto& face : allFaces) {
        if (!face.isBoundary()) continue;
        
        size_t P = face.ownerCell;  // Owner cell adjacent to boundary
        
        // Get boundary patch and condition
        const BoundaryPatch* patch = faceToPatchMap.at(face.id);
        const BoundaryData* bc = bcManager.getFieldBC(patch->patchName, "U");
        
        if (!bc) {
            // No BC specified for this patch - skip
            continue;
        }
        
        // Apply boundary condition based on type
        switch (bc->type) {
            case BCType::FIXED_VALUE:
                // For fixed value, directly set the cell value (can be improved with ghost cells)
                // This is a simple implementation - production code would use ghost cells
                U[P] = bc->getFixedVectorValue();
                break;
                
            case BCType::NO_SLIP:
                // No-slip wall: velocity = 0
                U[P] = Vector(0.0, 0.0, 0.0);
                break;
                
            case BCType::ZERO_GRADIENT:
                // Zero gradient: no modification needed as gradients are computed
                // The cell value remains as calculated by momentum equations
                break;
                
            case BCType::FIXED_GRADIENT:
                // Fixed gradient: modify based on gradient
                // Note: This is simplified - proper implementation would consider face distance
                break;
                
            default:
                break;
        }
    }
}

void SIMPLE::applyPressureBoundaryConditions() {
    // Build face-to-patch map for boundary lookup
    std::map<size_t, const BoundaryPatch*> faceToPatchMap;
    for (const auto& patch : bcManager.patches) {
        for (size_t i = patch.firstFaceIndex; i <= patch.lastFaceIndex; ++i) {
            faceToPatchMap[i] = &patch;
        }
    }
    
    // Apply pressure boundary conditions
    for (const auto& face : allFaces) {
        if (!face.isBoundary()) continue;
        
        size_t P = face.ownerCell;  // Owner cell adjacent to boundary
        
        // Get boundary patch and condition
        const BoundaryPatch* patch = faceToPatchMap.at(face.id);
        const BoundaryData* bc = bcManager.getFieldBC(patch->patchName, "p");
        
        if (!bc) {
            // No BC specified for this patch - skip
            continue;
        }
        
        // Apply boundary condition based on type
        switch (bc->type) {
            case BCType::FIXED_VALUE:
                // For fixed pressure outlet
                p[P] = bc->getFixedScalarValue();
                break;
                
            case BCType::ZERO_GRADIENT:
                // Zero gradient: no modification needed
                // The cell value remains as calculated by pressure correction
                break;
                
            case BCType::FIXED_GRADIENT:
                // Fixed gradient: typically used at inlets
                // Note: This is simplified - proper implementation would consider face distance
                break;
                
            default:
                break;
        }
    }
}

Scalar SIMPLE::calculateMassImbalance() const {
    Scalar totalImbalance = 0.0;
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar cellImbalance = 0.0;
        
        // Sum mass fluxes for each cell
        for (size_t j = 0; j < allCells[i].faceIndices.size(); ++j) {
            size_t faceIdx = allCells[i].faceIndices[j];
            int sign = allCells[i].faceSigns[j];
            cellImbalance += sign * massFlux[faceIdx];
        }
        
        totalImbalance += std::abs(cellImbalance);
    }
    
    return totalImbalance / allCells.size();
}

Scalar SIMPLE::calculateVelocityResidual() const {
    // Calculate RMS of velocity magnitude
    Scalar sumSq = 0.0;
    for (size_t i = 0; i < allCells.size(); ++i) {
        sumSq += U[i].magnitudeSquared();
    }
    return std::sqrt(sumSq / allCells.size());
}

Scalar SIMPLE::calculatePressureResidual() const {
    // Calculate RMS of pressure correction
    Scalar sumSq = 0.0;
    for (size_t i = 0; i < allCells.size(); ++i) {
        sumSq += p_prime[i] * p_prime[i];
    }
    return std::sqrt(sumSq / allCells.size());
}

// Setter methods
void SIMPLE::setRelaxationFactors(Scalar alpha_U_new, Scalar alpha_p_new) {
    alpha_U = alpha_U_new;
    alpha_p = alpha_p_new;
}

void SIMPLE::setConvergenceTolerance(Scalar tol) {
    tolerance = tol;
}

void SIMPLE::setMaxIterations(int maxIter) {
    maxIterations = maxIter;
}

void SIMPLE::enableTurbulenceModeling(bool enable) {
    enableTurbulence = enable;
    if (enable) {
        std::cout << "k-omega SST turbulence modeling enabled." << std::endl;
    } else {
        std::cout << "Laminar flow (turbulence modeling disabled)." << std::endl;
    }
}

// Turbulence getters
const ScalarField* SIMPLE::getTurbulentKineticEnergy() const {
    if (enableTurbulence && turbulenceModel) {
        return &(turbulenceModel->getK());
    }
    return nullptr;
}

const ScalarField* SIMPLE::getSpecificDissipationRate() const {
    if (enableTurbulence && turbulenceModel) {
        return &(turbulenceModel->getOmega());
    }
    return nullptr;
}

const ScalarField* SIMPLE::getTurbulentViscosity() const {
    if (enableTurbulence && turbulenceModel) {
        return &(turbulenceModel->getTurbulentViscosity());
    }
    return nullptr;
}

const ScalarField* SIMPLE::getWallDistance() const {
    if (enableTurbulence && turbulenceModel) {
        return &(turbulenceModel->getWallDistance());
    }
    return nullptr;
}

// Velocity component getters
ScalarField SIMPLE::getVelocityX() const {
    ScalarField U_x("U_x", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_x[i] = U[i].x;
    }
    return U_x;
}

ScalarField SIMPLE::getVelocityY() const {
    ScalarField U_y("U_y", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_y[i] = U[i].y;
    }
    return U_y;
}

ScalarField SIMPLE::getVelocityZ() const {
    ScalarField U_z("U_z", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_z[i] = U[i].z;
    }
    return U_z;
}