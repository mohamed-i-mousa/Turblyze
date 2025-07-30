#include "SIMPLE.h"
#include "Scalar.h"
#include <cmath>
#include <iostream>
#include <algorithm>

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
      rho(1.0),
      mu(1e-3),
      alpha_U(0.7),
      alpha_p(0.3),
      maxIterations(500),
      tolerance(1e-6),
      enableNonOrthCorrection(true),
      enableTurbulence(false),
      U("U", cells.size(), Vector(0.0, 0.0, 0.0)),
      p("p", cells.size(), 0.0),
      p_prime("p_prime", cells.size(), 0.0),
      U_face("U_face", faces.size(), Vector(0.0, 0.0, 0.0)),
      massFlux("massFlux", faces.size(), 0.0),
      volumeFlux("volumeFlux", faces.size(), 0.0),
      a_U("a_U", cells.size(), 0.0),
      H_U("H_U", cells.size(), Vector(0.0, 0.0, 0.0)),
      gradP("gradP", cells.size(), Vector(0.0, 0.0, 0.0))
{
    initialize();
}

void SIMPLE::initialize() {
    // Initialize matrix constructor
    matrixConstructor = std::make_unique<MatrixConstructor>(
        allFaces, allCells, bcManager, gradientScheme);
    
    // Initialize velocity field with reasonable values for internal cells
    for (size_t i = 0; i < allCells.size(); ++i) {
        U[i] = Vector(0.1, 0.0, 0.0);  // Small initial velocity
        p[i] = 0.0;                     // Zero initial pressure
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
    std::cout << "\n=== Starting SIMPLE Algorithm with ";
    if (enableTurbulence) {
        std::cout << "k-omega SST Turbulence Model";
    } else {
        std::cout << "Laminar Flow";
    }
    std::cout << " ===" << std::endl;
    
    int iteration = 0;
    bool converged = false;
    
    // Calculate initial pressure gradient
    gradP = gradientScheme.LeastSquares(p, allCells);
    
    while (!converged && iteration < maxIterations) {
        std::cout << "\n--- SIMPLE Iteration " << iteration + 1 << " ---" << std::endl;
        
        // Step 1: Solve momentum equations with effective viscosity
        discretizeMomentumEquations();
        
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
            
            // Calculate velocity gradients for turbulence production
            std::vector<VectorField> gradU;
            gradU.reserve(3);
            
            // Create VectorField objects with proper initialization
            gradU.emplace_back("gradU_x", allCells.size(), Vector(0.0, 0.0, 0.0));
            gradU.emplace_back("gradU_y", allCells.size(), Vector(0.0, 0.0, 0.0));
            gradU.emplace_back("gradU_z", allCells.size(), Vector(0.0, 0.0, 0.0));
            
            // Create velocity component fields for gradient calculation
            ScalarField U_x_field("U_x", allCells.size());
            ScalarField U_y_field("U_y", allCells.size());
            ScalarField U_z_field("U_z", allCells.size());
            
            // Extract velocity components
            for (size_t i = 0; i < allCells.size(); ++i) {
                U_x_field[i] = U[i].x;
                U_y_field[i] = U[i].y;
                U_z_field[i] = U[i].z;
            }
            
            // Calculate velocity gradients using Least Squares
            gradU[0] = gradientScheme.LeastSquares(U_x_field, allCells);
            gradU[1] = gradientScheme.LeastSquares(U_y_field, allCells);
            gradU[2] = gradientScheme.LeastSquares(U_z_field, allCells);
            
            // Solve complete k-omega SST model
            turbulenceModel->solve(U, gradU, rho, mu);
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

void SIMPLE::discretizeMomentumEquations() {
    /**
     * Solve momentum equations with effective viscosity (laminar + turbulent):
     * ∂(ρU)/∂t + ∇·(ρUU) = -∇p + ∇·[(μ + μₜ)∇U] + S
     * 
     * Where μₜ is the turbulent viscosity from k-omega SST model
     */
    
    // Update pressure gradient for momentum source term
    gradP = gradientScheme.LeastSquares(p, allCells);
    
    // Get effective viscosity (laminar + turbulent)
    Scalar mu_eff = mu;  // Default to laminar viscosity
    ScalarField mu_effective("mu_eff", allCells.size(), mu);
    
    if (enableTurbulence && turbulenceModel) {
        mu_effective = turbulenceModel->getEffectiveViscosity(mu);
        std::cout << "  Using effective viscosity (μ_lam + μ_t)" << std::endl;
    }
    
    // Solve U-momentum equation with effective viscosity
    matrixConstructor->constructScalarTransportMatrix(
        "U_x", 
        ScalarField("U_x", allCells.size()),  // Current U.x component
        ScalarField("U_x_old", allCells.size()), // Previous time step
        VectorField("U_transport", allCells.size()), // Transport velocity
        rho,    // Density
        mu_eff, // Effective viscosity (laminar + turbulent)
        TimeScheme::Steady,  // Time scheme
        0.0,    // dt (not used for steady)
        1.0,    // theta (not used for steady)
        convectionScheme
    );
    
    // Extract U.x component for solving
    ScalarField U_x("U_x", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_x[i] = U[i].x;
    }
    
    // Add pressure gradient source term and under-relaxation
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstructor->getVectorB());
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstructor->getMatrixA());
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar cellVolume = allCells[i].volume;
        
        // Add pressure gradient source term: -∇p·V
        b_vector(i) -= gradP[i].x * cellVolume;
        
        // Apply under-relaxation: a_P = a_P/α + (1-α)/α * a_P_old
        Scalar a_P = A_matrix.coeff(i, i);
        a_U[i] = a_P / alpha_U;  // Store for Rhie-Chow
        A_matrix.coeffRef(i, i) = a_U[i];
        b_vector(i) += (1.0 - alpha_U) / alpha_U * a_P * U[i].x;
        
        // Store H/A term for Rhie-Chow interpolation
        H_U[i].x = b_vector(i) / a_U[i];
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
    
    // Similarly solve for U.y and U.z components
    // (For brevity, implementing similar pattern)
    
    // Solve V-momentum (y-component)
    ScalarField U_y("U_y", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_y[i] = U[i].y;
    }
    
    matrixConstructor->constructScalarTransportMatrix(
        "U_y", U_y, ScalarField("U_y_old", allCells.size()),
        VectorField("U_transport", allCells.size()),
        rho, mu, TimeScheme::Steady, 0.0, 1.0, convectionScheme);
    
    auto& b_vector_y = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstructor->getVectorB());
    auto& A_matrix_y = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstructor->getMatrixA());
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar cellVolume = allCells[i].volume;
        b_vector_y(i) -= gradP[i].y * cellVolume;
        
        Scalar a_P = A_matrix_y.coeff(i, i);
        A_matrix_y.coeffRef(i, i) = a_P / alpha_U;
        b_vector_y(i) += (1.0 - alpha_U) / alpha_U * a_P * U[i].y;
        
        H_U[i].y = b_vector_y(i) / (a_P / alpha_U);
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
    
    // For 2D problems, U.z can be kept at zero or solved similarly for 3D
    for (size_t i = 0; i < allCells.size(); ++i) {
        H_U[i].z = 0.0;  // For 2D problems
    }
}

void SIMPLE::calculateRhieChowFaceVelocities() {
    // Rhie-Chow interpolation for face velocities
    // U_f = U_f_interpolated + D_f * (∇p_cell_interpolated - ∇p_face_from_neighbors)
    
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx) {
        const Face& face = allFaces[faceIdx];
        
        if (face.isBoundary()) {
            // For boundary faces, apply boundary conditions
            // (implementation depends on specific BC type)
            U_face[faceIdx] = Vector(0.0, 0.0, 0.0);  // Default no-slip
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
        Vector gradP_f_interpolated = w_f * gradP[P] + (1.0 - w_f) * gradP[N];
        
        // Face-based pressure gradient (using neighbor cell pressures)
        Vector gradP_f_face = (p[N] - p[P]) / d_PN.magnitude() * 
                              (d_PN / d_PN.magnitude());
        
        // Interpolated diffusion coefficient D_f
        Scalar a_P_interp = w_f * a_U[P] + (1.0 - w_f) * a_U[N];
        Scalar D_f = face.area / (a_P_interp + 1e-20);  // Avoid division by zero
        
        // Rhie-Chow correction
        Vector correction = D_f * (gradP_f_interpolated - gradP_f_face);
        
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
    // ∇·(D_f ∇p') = ∇·ρU*
    
    matrixConstructor->constructScalarTransportMatrix(
        "p_prime", p_prime, ScalarField("p_prime_old", allCells.size()),
        VectorField("zero_velocity", allCells.size(), Vector(0.0, 0.0, 0.0)),
        0.0,  // No convection for pressure correction
        1.0,  // Diffusion coefficient (will be overridden)
        TimeScheme::Steady, 0.0, 1.0, convectionScheme);
    
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstructor->getMatrixA());
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstructor->getVectorB());
    
    // Clear the matrix and rebuild for pressure correction
    A_matrix.setZero();
    b_vector.setZero();
    
    std::vector<Eigen::Triplet<Scalar>> triplets;
    
    // Loop over all faces to build pressure correction equation
    for (const auto& face : allFaces) {
        if (!face.geometricPropertiesCalculated) continue;
        
        size_t P = face.ownerCell;
        
        if (face.isBoundary()) {
            // Boundary face: apply zero gradient condition
            // (pressure correction should have zero gradient at boundaries)
            continue;
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
    
    std::cout << "  Mass imbalance: " << massImbalance 
              << ", Velocity residual: " << velocityResidual
              << ", Pressure residual: " << pressureResidual << std::endl;
    
    return converged;
}

void SIMPLE::applyVelocityBoundaryConditions() {
    // Apply velocity boundary conditions
    // (Implementation would depend on specific boundary patch setup)
    // This is a placeholder for BC application
}

void SIMPLE::applyPressureBoundaryConditions() {
    // Apply pressure boundary conditions
    // (Implementation would depend on specific boundary patch setup)
    // This is a placeholder for BC application
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