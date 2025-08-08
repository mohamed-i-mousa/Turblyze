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
               const ConvectionScheme& convScheme)
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
      U_prev("U_prev", cells.size(), Vector(0.0, 0.0, 0.0)),
      U_face_prev("U_face_prev", faces.size(), Vector(0.0, 0.0, 0.0)),
      a_U("a_U", cells.size(), 0.0),
      H_U("H_U", cells.size(), Vector(0.0, 0.0, 0.0)),
      gradP("gradP", cells.size(), Vector(0.0, 0.0, 0.0))
{
    initialize(Vector(0.0, 0.0, -0.1), 0.0);
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
    
    std::cout << "SIMPLE algorithm initialized with " << allCells.size() 
              << " cells and " << allFaces.size() << " faces." << std::endl;
}

void SIMPLE::solve() {

    std::cout << "\n=== Starting SIMPLE Loop" << " ===" << std::endl;
    
    int iteration = 0;
    bool converged = false;
    
    while (!converged && iteration < maxIterations) {
        std::cout << "\n--- SIMPLE Iteration " << iteration + 1 << " ---" << std::endl;

        U_prev = U;
        U_face_prev = U_face;
   
        // Refresh caches (gradients, mdot) for this iteration
        matrixConstruct->refreshIterationCaches(p, U, rho, turbulenceModel.get());
        gradP = matrixConstruct->gradP;
        
        solveMomentumEquations();
            
        calculateRhieChowFaceVelocities();
        calculateMassFluxes();
        
        solvePressureCorrection();
        
        correctVelocity();
        correctMassFluxes();
        correctPressure();
        
        if (enableTurbulence && turbulenceModel) {
            std::cout << "  Solving turbulence equations..." << std::endl;
            
            turbulenceModel->solve(U,
                                 matrixConstruct->gradUx,
                                 matrixConstruct->gradUy,
                                 matrixConstruct->gradUz,
                                 rho, mu);
        }
        
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

/*
 * Solve momentum equations:
 * ∂(ρU)/∂t + ∇·(ρUU) = -∇p + ∇·[(μ + μₜ)∇U] + S
  */
void SIMPLE::solveMomentumEquations() {

    // Reset diagonal accumulator
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_U[i] = 0.0;
    }
    
    Scalar mu_eff = mu;  // Laminar flow
    
    if (enableTurbulence && turbulenceModel) {
        // For turbulent flow, use average effective viscosity
        ScalarField mu_effective = turbulenceModel->getEffectiveViscosity(mu);
        Scalar sum_mu = 0.0;
        for (size_t i = 0; i < allCells.size(); ++i) {
            sum_mu += mu_effective[i];
        }
        mu_eff = sum_mu / allCells.size();
        std::cout << "  Using averaged effective viscosity (μ_lam + μ_t): " << mu_eff << std::endl;
    }
    
    // ********************************************************** //
    // *********************** x-momentum *********************** //
    // ********************************************************** //

    // Extract U.x component
    ScalarField U_x("U_x", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_x[i] = U[i].x;
    }
    
    // Calculate pressure gradient source term: Sx = -(∂p/∂x) * V
    ScalarField Ux_source("U_x_source", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        Ux_source[i] = -gradP[i].x * allCells[i].volume;
    }

    // Build x-momentum matrix
    matrixConstruct->buildMomentumMatrix(
        "U_x", 
        U_x,
        ScalarField("U_x_old", allCells.size()),
        Ux_source,
        rho,
        mu_eff,
        TimeScheme::Steady,
        0.0,
        1.0,
        matrixConstruct->gradUx,
        matrixConstruct->gradUx_f,
        convectionScheme
    );
    
    // Apply implicit under-relaxation before solving (Patankar's relaxation)
    matrixConstruct->relax(alpha_U, U_x);
    
    // Get references to the assembled matrix and vector
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstruct->getVectorB());
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstruct->getMatrixA());
    
    // Store the diagonal coefficient for Rhie-Chow interpolation (pre-relaxation)
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_U[i] = A_matrix.coeff(i, i);
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
    
    // ********************************************************** //
    // *********************** y-momentum *********************** //
    // ********************************************************** //

    // Extract U.y component
    ScalarField U_y("U_y", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_y[i] = U[i].y;
    }
    
    // Calculate pressure gradient source term: Sy = -(∂p/∂y) * V
    ScalarField Uy_source("U_y_source", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        Uy_source[i] = -gradP[i].y * allCells[i].volume;
    }

    // Build y-momentum matrix
    matrixConstruct->buildMomentumMatrix(
        "U_y", U_y, ScalarField("U_y_old", allCells.size()),
        Uy_source,
        rho, mu_eff, 
        TimeScheme::Steady,
        0.0, 
        1.0, 
        matrixConstruct->gradUy,
        matrixConstruct->gradUy_f,
        convectionScheme);
    
    // Get references to the assembled matrix and vector
    auto& b_vector_y = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstruct->getVectorB());
    auto& A_matrix_y = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstruct->getMatrixA());
        
    // Apply implicit under-relaxation before solving (Patankar's relaxation)
    matrixConstruct->relax(alpha_U, U_y);
    
    // Update blended diagonal coefficient for Rhie-Chow interpolation
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_U[i] = 0.5 * (a_U[i] + A_matrix_y.coeff(i, i));
    }
    
    // Solve for U.y
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
    
    // ********************************************************** //
    // *********************** z-momentum *********************** //
    // ********************************************************** //

    // Extract U.z component
    ScalarField U_z("U_z", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        U_z[i] = U[i].z;
    }
    
    // Build pressure source term: Sz = -(∂p/∂z) * V
    ScalarField Uz_source("U_z_source", allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        Uz_source[i] = -gradP[i].z * allCells[i].volume;
    }

    // Build z-momentum matrix
    matrixConstruct->buildMomentumMatrix(
        "U_z", U_z, ScalarField("U_z_old", allCells.size()),
        Uz_source,
        rho, mu_eff, 
        TimeScheme::Steady,
        0.0, 
        1.0, 
        matrixConstruct->gradUz,
        matrixConstruct->gradUz_f,
        convectionScheme);
    
    // Get references to the assembled matrix and vector
    auto& b_vector_z = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstruct->getVectorB());
    auto& A_matrix_z = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstruct->getMatrixA());
        
    // Apply implicit under-relaxation before solving (Patankar's relaxation)
    matrixConstruct->relax(alpha_U, U_z);
    
    // Update blended diagonal coefficient for Rhie-Chow interpolation
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_U[i] = (a_U[i] + A_matrix_z.coeff(i, i)) * 0.5;
    }
    
    // Solve for U.z
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
}

void SIMPLE::calculateRhieChowFaceVelocities() {
    // Rhie-Chow interpolation for face velocities
    // U_f = U_f_interpolated + D_f * (∇p_cell_interpolated - ∇p_face_from_neighbors)
    //          + (1-alpha_U) * (m_f_interpolated - m_f_from_neighbors)
    //          + DT_f (U_f_interpolated - U_f_from_neighbors)
    
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx) {
        const Face& face = allFaces[faceIdx];
        
        if (face.isBoundary()) {
            // Use centralized boundary condition handling
            Vector U_b = bcManager.calculateBoundaryFaceVectorValue(face, U, allCells, "U");
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
        Scalar total_dist = d_P + d_N + 1e-20;
        
        // Correct distance-weighted interpolation coefficients
        Scalar w_P = d_N / total_dist; // weight for owner cell P
        Scalar w_N = d_P / total_dist; // weight for neighbour cell N
        
        // Standard interpolated face velocity
        Vector U_f_interpolated = w_P * U[P] + w_N * U[N];
        
        // Interpolated pressure gradient
        Vector gradP_f_interpolated = w_P * gradP[P] + w_N * gradP[N];
        
        // Use the line-of-centres unit vector for the face pressure gradient
        Vector e_PN = d_PN / (d_PN.magnitude() + 1e-20);
        
        // Interpolated diffusion coefficient D_f (Rhie-Chow)
        Scalar a_P_interp = w_P * a_U[P] + w_N * a_U[N];
        Vector S_f = face.normal * face.area;
        Scalar dPN = d_PN.magnitude() + 1e-20;
        Scalar En = std::abs(dot(S_f, d_PN / dPN));
        // Use consistent coefficient including line-of-centres distance
        Scalar D_f = En / ((a_P_interp) * dPN + 1e-20);
        
        // Rhie-Chow correction: align with cell-centre line
        Vector gradP_f_face = ((p[N] - p[P]) / (dPN)) * e_PN;
        Vector correction = D_f * (gradP_f_interpolated - gradP_f_face);
        
        // Final face velocity with under-relaxation effect from previous iteration
        Vector U_f = U_f_interpolated - correction;

        if (hasPrevIterData && !face.isBoundary()) {
            // Compute previous-iteration interpolated face velocity using previous cell values
            Vector U_f_interpolated_prev = w_P * U_prev[P] + w_N * U_prev[N];
            // Add (1 - alpha_U) * (v_f_prev - v_f_avg_prev)
            U_f = U_f + (S(1.0) - alpha_U) * (U_face_prev[faceIdx] - U_f_interpolated_prev);
        }

        U_face[faceIdx] = U_f;
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

void SIMPLE::calculateFaceFluxes() {
    // Backward-compatible wrapper: compute face velocities then fluxes
    // Ensures that callers of calculateFaceFluxes() get up-to-date mass/volume fluxes
    calculateRhieChowFaceVelocities();
    calculateMassFluxes();
}

void SIMPLE::solvePressureCorrection() {
    // Build pressure correction equation
    // ∇·(D_f ∇p') = ∇·ρU  >> ∇·(D_f ∇p') = -∑ m*_f
    matrixConstruct->buildPressureMatrix(massFlux, a_U, rho);

    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(matrixConstruct->getMatrixA());
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(matrixConstruct->getVectorB());

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
    // Velocity correction: U = U* - (1/a_U) * ∇p'
    VectorField gradP_prime = gradientScheme.LeastSquares(p_prime, allCells);
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Vector correction = (S(1.0) / (a_U[i] + 1e-20)) * gradP_prime[i];
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

void SIMPLE::correctMassFluxes() {
    // Update mass flux on faces to be consistent with corrected pressure
    // m_f = m_f* - D_f * (p'_N - p'_P)
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx) {
        const Face& face = allFaces[faceIdx];
        if (face.isBoundary()) continue; // boundary handled by BCs

        size_t P = face.ownerCell;
        size_t N = face.neighbourCell.value();

        Vector d_PN = allCells[N].centroid - allCells[P].centroid;
        Vector d_Pf = face.centroid - allCells[P].centroid;
        Vector d_Nf = face.centroid - allCells[N].centroid;
        Scalar d_P = d_Pf.magnitude();
        Scalar d_N = d_Nf.magnitude();
        Scalar total_dist = d_P + d_N + 1e-20;
        Scalar w_P = d_N / total_dist;
        Scalar w_N = d_P / total_dist;
        Scalar a_face = w_P * a_U[P] + w_N * a_U[N];

        Vector S_f = face.normal * face.area;
        Vector e_PN = d_PN / (d_PN.magnitude() + 1e-20);
        Scalar E_mag = std::abs(dot(S_f, e_PN));
        Scalar dPN = d_PN.magnitude() + 1e-20;
        Scalar D_f = rho * (E_mag / (a_face + 1e-20)) / dPN;

        // Use alpha_p-scaled correction compatible with pressure update
        Scalar dp_prime = alpha_p * (p_prime[N] - p_prime[P]);
        massFlux[faceIdx] -= D_f * dp_prime;
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