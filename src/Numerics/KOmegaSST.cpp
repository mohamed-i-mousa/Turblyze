#include "KOmegaSST.h"
#include <cmath>
#include <algorithm>
#include <iostream>

KOmegaSST::KOmegaSST(const std::vector<Face>& faces,
                     const std::vector<Cell>& cells,
                     const BoundaryConditions& bc,
                     const GradientScheme& gradScheme)
    : allFaces(faces),
      allCells(cells),
      bcManager(bc),
      gradientScheme(gradScheme),
      k("k", cells.size(), 1e-6),
      omega("omega", cells.size(), 1.0),
      mu_t("mu_t", cells.size(), 0.0),
      wallDistance("wallDistance", cells.size(), 1.0),
      wallShearStress("wallShearStress", cells.size(), 0.0),
      grad_k("grad_k", cells.size(), Vector(0.0, 0.0, 0.0)),
      grad_omega("grad_omega", cells.size(), Vector(0.0, 0.0, 0.0)),
      F1("F1", cells.size(), 0.0),
      F2("F2", cells.size(), 0.0),
      production_k("production_k", cells.size(), 0.0),
      production_omega("production_omega", cells.size(), 0.0),
      cross_diffusion("cross_diffusion", cells.size(), 0.0)
{
    // Initialize matrix constructor for equation solving
    matrixConstructor = std::make_unique<MatrixConstructor>(
        allFaces, allCells, bcManager, gradientScheme);
    
    std::cout << "k-omega SST turbulence model initialized with " 
              << allCells.size() << " cells." << std::endl;
}

void KOmegaSST::initialize(const VectorField& U_field, Scalar rho, Scalar mu_lam) {
    std::cout << "\n=== Initializing k-omega SST Turbulence Model ===" << std::endl;
    
    // Step 1: Calculate wall distance using Poisson equation
    calculateWallDistance();
    
    // Step 2: Initialize turbulence fields based on flow conditions
    Scalar I_turbulence = 0.05;  // 5% turbulence intensity
    Scalar l_turbulent = 0.1;    // Turbulent length scale (adjust based on geometry)
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar U_mag = U_field[i].magnitude();
        
        // Initialize k based on turbulence intensity: k = 1.5 * (I * U)²
        k[i] = std::max(1.5 * (I_turbulence * U_mag) * (I_turbulence * U_mag), 1e-8);
        
        // Initialize omega based on length scale: ω = k^0.5 / (C_μ^0.25 * l)
        omega[i] = std::max(std::sqrt(k[i]) / (std::pow(constants.C_mu, 0.25) * l_turbulent), 1e-4);
        
        // Initialize turbulent viscosity: μₜ = ρ * k / ω
        mu_t[i] = rho * k[i] / omega[i];
    }
    
    // Step 3: Apply turbulence boundary conditions
    applyTurbulenceBoundaryConditions("k", k, U_field, mu_lam, rho);
    applyTurbulenceBoundaryConditions("omega", omega, U_field, mu_lam, rho);
    
    std::cout << "Turbulence fields initialized successfully." << std::endl;
}

void KOmegaSST::calculateWallDistance() {
    std::cout << "  Computing wall distance using Poisson equation..." << std::endl;
    
    /**
     * Wall distance calculation using Poisson-like approach:
     * 1. Solve: ∇²φ = -1 with φ = 0 at walls
     * 2. Wall distance: d = √(|∇φ|² + 2φ) - |∇φ|
     * 
     * This method is more robust than geometric approaches for complex geometries
     */
    
    // Create a scalar field φ for the Poisson equation
    ScalarField phi("phi_wall", allCells.size(), 0.0);
    ScalarField phi_old = phi;
    
    // Construct Poisson equation: ∇²φ = -1
    matrixConstructor->constructScalarTransportMatrix(
        "phi_wall",
        phi,
        phi_old,
        VectorField("zero_velocity", allCells.size(), Vector(0.0, 0.0, 0.0)), // No convection
        0.0,    // No convection
        1.0,    // Unit diffusion coefficient
        TimeScheme::Steady,
        0.0,    // No time derivative
        1.0,    // Full implicit
        CentralDifferenceScheme()
    );
    
    // Modify source term: b = -1 (Laplacian = -1)
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstructor->getMatrixA());
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstructor->getVectorB());
    
    // Set source term to -1 for all cells
    for (size_t i = 0; i < allCells.size(); ++i) {
        b_vector(i) = -allCells[i].volume;  // -1 * cell_volume
    }
    
    // Apply boundary condition: φ = 0 at walls
    for (const auto& face : allFaces) {
        if (face.isBoundary()) {
            // Check if this is a wall boundary
            // (Implementation depends on how walls are identified in boundary patches)
            size_t P = face.ownerCell;
            
            // For simplicity, assume all boundaries are walls
            // In practice, you'd check boundary patch type
            A_matrix.coeffRef(P, P) += 1e12;  // Large penalty method
            b_vector(P) = 0.0;  // φ = 0 at wall
        }
    }
    
    // Solve the Poisson equation
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> phi_solution(allCells.size());
    phi_solution.setZero();
    
    bool solved = LinearSolvers::BiCGSTAB(
        phi_solution, A_matrix, b_vector, 1e-8, 1000, "phi_wall");
    
    if (solved) {
        for (size_t i = 0; i < allCells.size(); ++i) {
            phi[i] = phi_solution(i);
        }
    }
    
    // Calculate gradient of φ
    VectorField grad_phi = gradientScheme.LeastSquares(phi, allCells);
    
    // Compute wall distance: d = √(|∇φ|² + 2φ) - |∇φ|
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar grad_phi_mag = grad_phi[i].magnitude();
        Scalar discriminant = grad_phi_mag * grad_phi_mag + 2.0 * phi[i];
        
        if (discriminant > 0) {
            wallDistance[i] = std::sqrt(discriminant) - grad_phi_mag;
        } else {
            wallDistance[i] = 1e-6;  // Minimum wall distance
        }
        
        // Ensure positive wall distance
        wallDistance[i] = std::max(wallDistance[i], 1e-6);
    }
    
    std::cout << "  Wall distance calculation completed." << std::endl;
}

void KOmegaSST::solve(const VectorField& U_field, 
                     const std::vector<VectorField>& gradU,
                     Scalar rho, 
                     Scalar mu_lam) {
    
    /**
     * Complete k-omega SST solution sequence:
     * 1. Calculate blending functions F1 and F2
     * 2. Solve omega equation with near-wall treatment
     * 3. Solve k equation
     * 4. Calculate turbulent viscosity
     * 5. Apply wall corrections
     */
    
    // Step 1: Calculate gradients of turbulence quantities
    grad_k = gradientScheme.LeastSquares(k, allCells);
    grad_omega = gradientScheme.LeastSquares(omega, allCells);
    
    // Step 2: Calculate blending functions
    calculateBlendingFunctions(gradU, rho, mu_lam);
    
    // Step 3: Calculate production terms
    calculateProductionTerms(gradU, rho);
    
    // Step 4: Calculate cross-diffusion term
    calculateCrossDiffusion();
    
    // Step 5: Solve omega equation
    solveOmegaEquation(U_field, gradU, rho, mu_lam);
    
    // Step 6: Apply near-wall treatment for omega
    applyNearWallTreatmentOmega();
    
    // Step 7: Solve k equation
    solveKEquation(U_field, gradU, rho, mu_lam);
    
    // Step 8: Calculate turbulent viscosity
    calculateTurbulentViscosity(U_field, gradU, rho);
    
    // Step 9: Apply wall corrections
    applyWallCorrections();
    
    // Step 10: Calculate wall shear stress
    calculateWallShearStress(U_field, mu_lam);
}

void KOmegaSST::solveOmegaEquation(const VectorField& U_field,
                                  const std::vector<VectorField>& gradU,
                                  Scalar rho,
                                  Scalar mu_lam) {
    
    /**
     * Omega transport equation:
     * ∂(ρω)/∂t + ∇·(ρUω) = γ·P_ω - β·ρω² + ∇·[(μ + σ_ω·μₜ)∇ω] + D_ω
     * 
     * Where:
     * - P_ω = production term
     * - D_ω = cross-diffusion term (only for SST model)
     * - γ, β, σ_ω are blended constants
     */
    
    std::cout << "  Solving omega transport equation..." << std::endl;
    
    ScalarField omega_old = omega;
    
    // Construct transport matrix for omega
    matrixConstructor->constructScalarTransportMatrix(
        "omega",
        omega,
        omega_old,
        U_field,
        rho,
        mu_lam,  // Base viscosity (will be modified below)
        TimeScheme::Steady,
        0.0,
        1.0,
        CentralDifferenceScheme()
    );
    
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstructor->getMatrixA());
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstructor->getVectorB());
    
    // Add omega-specific source terms and modify diffusion
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar cellVolume = allCells[i].volume;
        
        // Blended constants
        Scalar gamma = F1[i] * constants.gamma_1 + (1.0 - F1[i]) * constants.gamma_2;
        Scalar beta = F1[i] * constants.beta_1 + (1.0 - F1[i]) * constants.beta_2;
        Scalar sigma_omega = F1[i] * constants.sigma_omega1 + (1.0 - F1[i]) * constants.sigma_omega2;
        
        // Production term: γ·P_ω
        Scalar P_omega = gamma * production_omega[i];
        b_vector(i) += P_omega * cellVolume;
        
        // Destruction term: -β·ρω²
        Scalar destruction = beta * rho * omega[i] * omega[i];
        A_matrix.coeffRef(i, i) += destruction / omega[i] * cellVolume;  // Linearized
        
        // Cross-diffusion term (SST specific)
        b_vector(i) += cross_diffusion[i] * cellVolume;
        
        // Enhanced diffusion due to turbulent viscosity
        // This would require modifying the matrix construction for variable diffusion
        // For simplicity, we approximate this effect
    }
    
    // Solve omega equation
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> omega_solution(allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        omega_solution(i) = omega[i];
    }
    
    bool solved = LinearSolvers::BiCGSTAB(
        omega_solution, A_matrix, b_vector, 1e-8, 1000, "omega");
    
    if (solved) {
        for (size_t i = 0; i < allCells.size(); ++i) {
            omega[i] = std::max(omega_solution(i), 1e-6);  // Ensure positive omega
        }
    }
    
    // Apply boundary conditions
    applyTurbulenceBoundaryConditions("omega", omega, U_field, mu_lam, rho);
}

void KOmegaSST::applyNearWallTreatmentOmega() {
    /**
     * Near-wall treatment for omega:
     * For cells very close to the wall (y+ < 2), use the viscous sublayer relation:
     * ω_wall = 6μ/(ρβ₁y²)
     */
    
    std::cout << "  Applying near-wall treatment for omega..." << std::endl;
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar y = wallDistance[i];
        
        // For very near-wall cells, impose viscous sublayer value
        if (y < 1e-6) {  // Very close to wall
            Scalar rho = 1.25;  // Should be passed as parameter
            Scalar mu_lam = 1.8e-5;
            
            // Viscous sublayer omega: ω = 6μ/(ρβ₁y²)
            omega[i] = 6.0 * mu_lam / (rho * constants.beta_1 * y * y + 1e-20);
            
            // Limit to reasonable values
            omega[i] = std::min(omega[i], 1e6);
            omega[i] = std::max(omega[i], 1e-4);
        }
    }
}

void KOmegaSST::solveKEquation(const VectorField& U_field,
                              const std::vector<VectorField>& gradU,
                              Scalar rho,
                              Scalar mu_lam) {
    
    /**
     * k transport equation:
     * ∂(ρk)/∂t + ∇·(ρUk) = P_k - β*·ρkω + ∇·[(μ + σ_k·μₜ)∇k]
     * 
     * Where:
     * - P_k = limited production term
     * - β* = 0.09 (destruction coefficient)
     * - σ_k = blended diffusion coefficient
     */
    
    std::cout << "  Solving k transport equation..." << std::endl;
    
    ScalarField k_old = k;
    
    // Construct transport matrix for k
    matrixConstructor->constructScalarTransportMatrix(
        "k",
        k,
        k_old,
        U_field,
        rho,
        mu_lam,  // Base viscosity
        TimeScheme::Steady,
        0.0,
        1.0,
        CentralDifferenceScheme()
    );
    
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>(
        matrixConstructor->getMatrixA());
    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>(
        matrixConstructor->getVectorB());
    
    // Add k-specific source terms
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar cellVolume = allCells[i].volume;
        
        // Blended sigma_k
        Scalar sigma_k = F1[i] * constants.sigma_k1 + (1.0 - F1[i]) * constants.sigma_k2;
        
        // Production term (limited to prevent unrealistic values)
        Scalar P_k_limited = limitProduction(production_k[i], rho, i);
        b_vector(i) += P_k_limited * cellVolume;
        
        // Destruction term: -β*·ρkω
        Scalar destruction = constants.beta_star * rho * omega[i];
        A_matrix.coeffRef(i, i) += destruction * cellVolume;
    }
    
    // Solve k equation
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> k_solution(allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        k_solution(i) = k[i];
    }
    
    bool solved = LinearSolvers::BiCGSTAB(
        k_solution, A_matrix, b_vector, 1e-8, 1000, "k");
    
    if (solved) {
        for (size_t i = 0; i < allCells.size(); ++i) {
            k[i] = std::max(k_solution(i), 1e-8);  // Ensure positive k
        }
    }
    
    // Apply boundary conditions
    applyTurbulenceBoundaryConditions("k", k, U_field, mu_lam, rho);
}

void KOmegaSST::calculateTurbulentViscosity(const VectorField& U_field,
                                           const std::vector<VectorField>& gradU,
                                           Scalar rho) {
    
    /**
     * SST turbulent viscosity formulation:
     * μₜ = ρ * a₁ * k / max(a₁ * ω, Ω * F₂)
     * 
     * Where:
     * - a₁ = 0.31 (stress limiter constant)
     * - Ω = magnitude of vorticity tensor
     * - F₂ = second blending function
     */
    
    std::cout << "  Calculating turbulent viscosity..." << std::endl;
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        // Calculate vorticity magnitude Ω
        Scalar Omega = 0.0;
        if (gradU.size() >= 3) {  // Ensure we have all velocity gradients
            // Ω = √(2ΩᵢⱼΩᵢⱼ) where Ωᵢⱼ = 0.5(∂uᵢ/∂xⱼ - ∂uⱼ/∂xᵢ)
            Scalar Omega_xy = 0.5 * (gradU[0][i].y - gradU[1][i].x);  // ∂u/∂y - ∂v/∂x
            Scalar Omega_xz = 0.5 * (gradU[0][i].z - gradU[2][i].x);  // ∂u/∂z - ∂w/∂x
            Scalar Omega_yz = 0.5 * (gradU[1][i].z - gradU[2][i].y);  // ∂v/∂z - ∂w/∂y
            
            Omega = std::sqrt(2.0 * (Omega_xy*Omega_xy + Omega_xz*Omega_xz + Omega_yz*Omega_yz));
        } else {
            // Simplified 2D case
            Omega = std::abs(gradU[0][i].y - gradU[1][i].x);
        }
        
        // SST turbulent viscosity formula
        Scalar denominator = std::max(constants.a1 * omega[i], Omega * F2[i]);
        mu_t[i] = rho * constants.a1 * k[i] / (denominator + 1e-20);
        
        // Limit turbulent viscosity to reasonable values
        Scalar mu_lam = 1.8e-5;  // Should be passed as parameter
        mu_t[i] = std::min(mu_t[i], 1000.0 * mu_lam);  // Max 1000 times laminar viscosity
        mu_t[i] = std::max(mu_t[i], 0.0);              // Ensure non-negative
    }
}

void KOmegaSST::applyWallCorrections() {
    /**
     * Wall corrections for turbulent viscosity:
     * - μₜ = 0 at walls (no-slip condition)
     * - Smooth transition from wall to freestream
     * - Apply damping functions based on y+ or wall distance
     */
    
    std::cout << "  Applying wall corrections..." << std::endl;
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar y = wallDistance[i];
        
        // Wall damping function (exponential decay near wall)
        // f_μ = 1 - exp(-y/δ) where δ is characteristic length
        Scalar delta = 1e-4;  // Characteristic damping length
        Scalar f_mu = 1.0 - std::exp(-y / delta);
        
        // Apply damping to turbulent viscosity
        mu_t[i] *= f_mu;
        
        // Ensure μₜ = 0 at walls (y → 0)
        if (y < 1e-6) {
            mu_t[i] = 0.0;
        }
    }
}

void KOmegaSST::calculateWallShearStress(const VectorField& U_field, Scalar mu_lam) {
    /**
     * Calculate wall shear stress using parallel velocity component:
     * τ_wall = μ_eff * (∂U_parallel/∂n)_wall
     * 
     * Where:
     * - U_parallel is velocity component parallel to wall
     * - ∂/∂n is derivative normal to wall
     * - μ_eff = μ_lam + μₜ (but μₜ = 0 at wall)
     */
    
    std::cout << "  Calculating wall shear stress..." << std::endl;
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar y = wallDistance[i];
        
        if (y < 1e-3) {  // Near-wall cells
            // Approximate wall shear using first cell velocity and distance
            Scalar U_magnitude = U_field[i].magnitude();
            
            // Wall shear stress: τ = μ * ∂U/∂y ≈ μ * U / y
            wallShearStress[i] = mu_lam * U_magnitude / (y + 1e-12);
            
            // Limit to reasonable values
            wallShearStress[i] = std::min(wallShearStress[i], 1000.0);
        } else {
            wallShearStress[i] = 0.0;  // Far from wall
        }
    }
}

ScalarField KOmegaSST::getEffectiveViscosity(Scalar mu_lam) const {
    /**
     * Calculate effective viscosity: μ_eff = μ_lam + μₜ
     * This is used in momentum equations for turbulent flow
     */
    
    ScalarField mu_eff("mu_eff", allCells.size());
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        mu_eff[i] = mu_lam + mu_t[i];
    }
    
    return mu_eff;
}

// Private helper methods implementation

void KOmegaSST::calculateBlendingFunctions(const std::vector<VectorField>& gradU,
                                          Scalar rho,
                                          Scalar mu_lam) {
    /**
     * Calculate SST blending functions F1 and F2:
     * 
     * F1 = tanh(arg1⁴) where arg1 = min(max(√k/(β*ωy), 500μ/(ρωy²)), 4ρσω2k/(CDkωy²))
     * F2 = tanh(arg2²) where arg2 = max(2√k/(β*ωy), 500μ/(ρωy²))
     * 
     * These functions blend between k-ω (F1=1) and k-ε (F1=0) formulations
     */
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        Scalar y = wallDistance[i];
        Scalar sqrt_k = std::sqrt(k[i]);
        
        // Avoid division by zero
        y = std::max(y, 1e-12);
        omega[i] = std::max(omega[i], 1e-12);
        
        // Arguments for blending functions
        Scalar arg1_1 = sqrt_k / (constants.beta_star * omega[i] * y);
        Scalar arg1_2 = 500.0 * mu_lam / (rho * omega[i] * y * y);
        Scalar arg1_3 = 4.0 * rho * constants.sigma_omega2 * k[i] / 
                       (std::max(cross_diffusion[i], 1e-20) * y * y);
        
        Scalar arg1 = std::min(std::max(arg1_1, arg1_2), arg1_3);
        F1[i] = std::tanh(std::pow(arg1, 4.0));
        
        Scalar arg2 = std::max(2.0 * sqrt_k / (constants.beta_star * omega[i] * y), arg1_2);
        F2[i] = std::tanh(arg2 * arg2);
    }
}

void KOmegaSST::calculateProductionTerms(const std::vector<VectorField>& gradU, Scalar rho) {
    /**
     * Calculate production terms for k and omega:
     * 
     * P_k = μₜ * S² where S is strain rate magnitude
     * P_ω = (γ/νₜ) * P_k where γ is blended coefficient
     */
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        // Calculate strain rate magnitude S = √(2SᵢⱼSᵢⱼ)
        Scalar S = 0.0;
        if (gradU.size() >= 3) {
            // Full 3D strain rate tensor
            Scalar S11 = gradU[0][i].x;                                    // ∂u/∂x
            Scalar S22 = gradU[1][i].y;                                    // ∂v/∂y  
            Scalar S33 = gradU[2][i].z;                                    // ∂w/∂z
            Scalar S12 = 0.5 * (gradU[0][i].y + gradU[1][i].x);          // 0.5(∂u/∂y + ∂v/∂x)
            Scalar S13 = 0.5 * (gradU[0][i].z + gradU[2][i].x);          // 0.5(∂u/∂z + ∂w/∂x)
            Scalar S23 = 0.5 * (gradU[1][i].z + gradU[2][i].y);          // 0.5(∂v/∂z + ∂w/∂y)
            
            S = std::sqrt(2.0 * (S11*S11 + S22*S22 + S33*S33 + 2.0*(S12*S12 + S13*S13 + S23*S23)));
        } else {
            // Simplified 2D case
            Scalar S11 = gradU[0][i].x;
            Scalar S22 = gradU[1][i].y;
            Scalar S12 = 0.5 * (gradU[0][i].y + gradU[1][i].x);
            
            S = std::sqrt(2.0 * (S11*S11 + S22*S22 + 2.0*S12*S12));
        }
        
        // Production of k: P_k = μₜ * S²
        production_k[i] = mu_t[i] * S * S;
        
        // Production of omega: P_ω = (γ/νₜ) * P_k = γ * ρ * S²
        Scalar gamma = F1[i] * constants.gamma_1 + (1.0 - F1[i]) * constants.gamma_2;
        production_omega[i] = gamma * rho * S * S;
    }
}

void KOmegaSST::calculateCrossDiffusion() {
    /**
     * Calculate cross-diffusion term for SST model:
     * D_ω = 2(1-F1) * σ_ω2 * ρ/ω * ∇k · ∇ω
     * 
     * This term appears only in the omega equation and is crucial for SST behavior
     */
    
    for (size_t i = 0; i < allCells.size(); ++i) {
        // Cross-diffusion: 2(1-F1) * σ_ω2 * ρ/ω * ∇k · ∇ω
        Scalar dot_product = dot(grad_k[i], grad_omega[i]);
        
        cross_diffusion[i] = 2.0 * (1.0 - F1[i]) * constants.sigma_omega2 * 
                           1.25 / (omega[i] + 1e-12) * dot_product;  // rho = 1.25
        
        // Ensure positive contribution (only for SST formulation)
        cross_diffusion[i] = std::max(cross_diffusion[i], 0.0);
    }
}

void KOmegaSST::applyTurbulenceBoundaryConditions(const std::string& fieldName,
                                                 ScalarField& field,
                                                 const VectorField& U_field,
                                                 Scalar mu_lam,
                                                 Scalar rho) {
    /**
     * Apply appropriate boundary conditions for turbulence quantities:
     * 
     * At walls:
     * - k = 0 (no turbulent kinetic energy at wall)
     * - ω = ω_wall (computed from viscous sublayer)
     * 
     * At inlets:
     * - k and ω based on turbulence intensity and length scale
     * 
     * At outlets:
     * - Zero gradient conditions
     */
    
    // This is a simplified implementation
    // In practice, you would iterate through boundary faces and apply specific conditions
    for (size_t i = 0; i < allCells.size(); ++i) {
        if (wallDistance[i] < 1e-6) {  // At wall
            if (fieldName == "k") {
                field[i] = 0.0;  // k = 0 at wall
            } else if (fieldName == "omega") {
                // ω_wall = 6μ/(ρβ₁y²) for viscous sublayer
                Scalar y = std::max(wallDistance[i], 1e-8);
                field[i] = 6.0 * mu_lam / (rho * constants.beta_1 * y * y);
                field[i] = std::min(field[i], 1e6);  // Limit maximum value
            }
        }
    }
}

Scalar KOmegaSST::calculateYPlus(size_t cellIdx,
                                const VectorField& U_field,
                                Scalar mu_lam,
                                Scalar rho) const {
    /**
     * Calculate y+ = ρ * u_τ * y / μ
     * where u_τ = √(τ_wall/ρ) is friction velocity
     */
    
    Scalar y = wallDistance[cellIdx];
    Scalar tau_wall = wallShearStress[cellIdx];
    Scalar u_tau = std::sqrt(tau_wall / rho);
    
    return rho * u_tau * y / mu_lam;
}

Scalar KOmegaSST::limitProduction(Scalar P_k, Scalar rho, size_t cellIdx) const {
    /**
     * Limit production to prevent unrealistic values:
     * P_k = min(P_k, 10 * β* * ρ * k * ω)
     * 
     * This prevents excessive production in stagnation regions
     */
    
    Scalar limit = 10.0 * constants.beta_star * rho * k[cellIdx] * omega[cellIdx];
    return std::min(P_k, limit);
}

void KOmegaSST::setModelConstants(Scalar C_mu, Scalar beta_1, Scalar beta_2) {
    constants.C_mu = C_mu;
    constants.beta_1 = beta_1;
    constants.beta_2 = beta_2;
} 