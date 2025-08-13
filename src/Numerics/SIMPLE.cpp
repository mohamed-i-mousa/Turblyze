#include "SIMPLE.h"
#include "Scalar.h"
#include "massFlowRate.h"
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
      a_Ux("a_Ux", cells.size(), 0.0),
      a_Uy("a_Uy", cells.size(), 0.0),
      a_Uz("a_Uz", cells.size(), 0.0),
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

    // Initialize face velocities
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx) {
        const Face& face = allFaces[faceIdx];
        if (face.isBoundary()) {
            U_face[faceIdx] = bcManager.calculateBoundaryFaceVectorValue(face, U, "U");
        } else {
            U_face[faceIdx] = linearInterpolation(face, U);
        }
        U_face_prev[faceIdx] = U_face[faceIdx];
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
        
        // Use previous-iteration face velocities to assemble momentum with consistent mdot
        // Assemble momentum with mdot from previous iteration to avoid coupling lag
        {
            // Build mdot from previous iteration's face velocities for consistency
            FaceFluxField mdot_prev("mdot_prev", allFaces.size(), 0.0);
            for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx) {
                const Face& face = allFaces[faceIdx];
                const Vector S_f = face.normal * face.area;
                mdot_prev[faceIdx] = rho * dot(U_face_prev[faceIdx], S_f);
            }
            matrixConstruct->setFaceMassFluxes(mdot_prev);
        }
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

    // Reset diagonals accumulators
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_Ux[i] = 0.0;
        a_Uy[i] = 0.0;
        a_Uz[i] = 0.0;
    }
    
    // Build per-cell effective viscosity μ_eff = μ + μ_t (spatially varying when turbulence enabled)
    ScalarField mu_eff("mu_eff", allCells.size());
    if (enableTurbulence && turbulenceModel) {
        mu_eff = turbulenceModel->getEffectiveViscosity(mu);
    } else {
        for (size_t i = 0; i < allCells.size(); ++i) mu_eff[i] = mu;
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
    
    // Get references to the assembled matrix and vector
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstruct->getVectorB());
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstruct->getMatrixA());

    // Store the diagonal coefficient for Rhie-Chow interpolation
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_Ux[i] = A_matrix.coeff(i, i);
    }

    // Apply implicit under-relaxation before solving (Patankar's relaxation)
    matrixConstruct->relax(alpha_U, U_x);
    
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

    // Store y-diagonal
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_Uy[i] = A_matrix_y.coeff(i, i);
    }

    // Apply implicit under-relaxation before solving (Patankar's relaxation)
    matrixConstruct->relax(alpha_U, U_y);
    
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

    // Store z-diagonal
    for (size_t i = 0; i < allCells.size(); ++i) {
        a_Uz[i] = A_matrix_z.coeff(i, i);
    }

    // Apply implicit under-relaxation before solving (Patankar's relaxation)
    matrixConstruct->relax(alpha_U, U_z);
    
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
    // U_f = U_f_linearInterpolated
    //       + D_f * (∇p_face_linearInterpolated - ∇p_face_from_cache)
    //       + (1 - alpha_U) * (U_f_previousIter - U_f_previousIter_linearInterpolated)

    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx) {
        const Face& face = allFaces[faceIdx];

        if (face.isBoundary()) {
            // Centralized boundary condition handling for velocity
            U_face[faceIdx] = bcManager.calculateBoundaryFaceVectorValue(face, U, "U");
            continue;
        }

        const size_t P = face.ownerCell;
        const size_t N = face.neighbourCell.value();

        // Linear interpolation weights
        Scalar w_P = 0.0, w_N = 0.0;
        computeLinearWeights(face, w_P, w_N);

        // Linear-interpolated values
        const Vector U_f_lin = w_P * U[P] + w_N * U[N];
        const Vector gradP_f_lin = w_P * gradP[P] + w_N * gradP[N];

        // Cached face gradient from gradient scheme (more accurate face value)
        const Vector gradP_f_cache = matrixConstruct->gradP_f[faceIdx];

        // Interpolated momentum diagonals on face per component
        const Scalar a_face_x = w_P * a_Ux[P] + w_N * a_Ux[N];
        const Scalar a_face_y = w_P * a_Uy[P] + w_N * a_Uy[N];
        const Scalar a_face_z = w_P * a_Uz[P] + w_N * a_Uz[N];

        // Minimum-correction D_f with non-orthogonal consideration (face-velocity RC, no rho)
        const Vector d_PN = allCells[N].centroid - allCells[P].centroid;
        const Vector S_f = face.normal * face.area;
        const Scalar denom = d_PN.magnitudeSquared() + 1e-20;
        const Scalar numer = d_PN.x * S_f.x / (a_face_x + 1e-20)
                           + d_PN.y * S_f.y / (a_face_y + 1e-20)
                           + d_PN.z * S_f.z / (a_face_z + 1e-20);
        const Scalar D_f = numer / denom;

        // Core Rhie-Chow correction using face gradient cache
        Vector U_f = U_f_lin + D_f * (gradP_f_lin - gradP_f_cache);

        // Previous-iteration face under-relaxation term
        const Vector U_f_lin_prev = w_P * U_prev[P] + w_N * U_prev[N];
        U_f += (S(1.0) - alpha_U) * (U_face_prev[faceIdx] - U_f_lin_prev);

        U_face[faceIdx] = U_f;
    }
}

// ----- Linear interpolation helpers ----- //
void SIMPLE::computeLinearWeights(const Face& face, Scalar& w_P, Scalar& w_N) const {
    if (face.isBoundary()) {
        w_P = S(1.0);
        w_N = 0.0;
        return;
    }
    const Scalar d_P = face.d_Pf_mag;
    const Scalar d_N = face.d_Nf_mag.value();
    const Scalar total = d_P + d_N + 1e-20;
    w_P = d_N / total;
    w_N = d_P / total;
}

Vector SIMPLE::linearInterpolation(const Face& face, const VectorField& cellField) const {
    if (face.isBoundary()) {
        return cellField[face.ownerCell];
    }
    Scalar w_P = 0.0, w_N = 0.0;
    computeLinearWeights(face, w_P, w_N);
    const size_t P = face.ownerCell;
    const size_t N = face.neighbourCell.value();
    return w_P * cellField[P] + w_N * cellField[N];
}

Scalar SIMPLE::linearInterpolation(const Face& face, const ScalarField& cellField) const {
    if (face.isBoundary()) {
        return cellField[face.ownerCell];
    }
    Scalar w_P = 0.0, w_N = 0.0;
    computeLinearWeights(face, w_P, w_N);
    const size_t P = face.ownerCell;
    const size_t N = face.neighbourCell.value();
    return w_P * cellField[P] + w_N * cellField[N];
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
    matrixConstruct->buildPressureMatrix(massFlux, a_Ux, a_Uy, a_Uz, rho);

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
        // Use a representative diagonal; choose the maximum for stability
        Scalar a_rep = std::max(a_Ux[i], std::max(a_Uy[i], a_Uz[i]));
        Vector correction = (S(1.0) / (a_rep + 1e-20)) * gradP_prime[i];
        U[i] = U[i] - correction;
    }
}

void SIMPLE::correctPressure() {
    // Pressure correction: p = p + α_p * p'
    // Track p' RMS before reset to report meaningful pressure residual
    Scalar sumSq = 0.0;
    for (size_t i = 0; i < allCells.size(); ++i) {
        sumSq += p_prime[i] * p_prime[i];
    }
    lastPressureCorrectionRMS = std::sqrt(sumSq / allCells.size());

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
        Scalar d_Pf = face.d_Pf_mag;
        Scalar d_Nf = face.d_Nf_mag.value();
        Scalar total_dist = d_Pf + d_Nf + 1e-20;
        Scalar w_P = d_Nf / total_dist;
        Scalar w_N = d_Pf / total_dist;
        Scalar a_face_x = w_P * a_Ux[P] + w_N * a_Ux[N];
        Scalar a_face_y = w_P * a_Uy[P] + w_N * a_Uy[N];
        Scalar a_face_z = w_P * a_Uz[P] + w_N * a_Uz[N];

        Vector S_f = face.normal * face.area;
        Scalar denom = d_PN.magnitudeSquared() + 1e-20;
        Scalar numer = d_PN.x * S_f.x / (a_face_x + 1e-20)
                     + d_PN.y * S_f.y / (a_face_y + 1e-20)
                     + d_PN.z * S_f.z / (a_face_z + 1e-20);
        Scalar D_f = rho * (numer / denom);

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
    // Dimensionless normalized continuity residual per cell, averaged
    Scalar totalNormImbalance = 0.0;
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar net = 0.0;
        Scalar sumAbs = 0.0;
        for (size_t j = 0; j < allCells[i].faceIndices.size(); ++j) {
            const size_t faceIdx = allCells[i].faceIndices[j];
            const int sign = allCells[i].faceSigns[j];
            const Scalar mf = massFlux[faceIdx];
            net += sign * mf;
            sumAbs += std::abs(mf);
        }
        const Scalar denom = sumAbs + 1e-30;
        totalNormImbalance += std::abs(net) / denom;
    }
    return totalNormImbalance / std::max<size_t>(1, allCells.size());
}

Scalar SIMPLE::calculateVelocityResidual() const {
    // Normalized delta-based residual: ||U - U_prev||_2 / (||U_prev||_2 + eps)
    Scalar num = 0.0;
    Scalar den = 0.0;
    for (size_t i = 0; i < allCells.size(); ++i) {
        const Vector d = U[i] - U_prev[i];
        num += d.magnitudeSquared();
        den += U_prev[i].magnitudeSquared();
    }
    num = std::sqrt(num + 1e-30);
    den = std::sqrt(den + 1e-30);
    return num / den;
}

Scalar SIMPLE::calculatePressureResidual() const {
    // Normalize p' RMS by RMS(p)
    Scalar sumP2 = 0.0;
    for (size_t i = 0; i < allCells.size(); ++i) sumP2 += p[i] * p[i];
    const Scalar p_rms = std::sqrt(sumP2 / (allCells.size() + 1e-30)) + 1e-30;
    return lastPressureCorrectionRMS / p_rms;
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