/******************************************************************************
 * @file KOmegaSST.cpp
 * @brief Implementation of k-omega SST turbulence model
 *****************************************************************************/

#include "KOmegaSST.h"
#include "Matrix.h"
#include "ConvectionScheme.h"
#include "GradientScheme.h"
#include <cmath>
#include <algorithm>
#include <iostream>

KOmegaSST::KOmegaSST
(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const BoundaryConditions& bc,
    const GradientScheme& gradientScheme
)
    : k("k", cells.size(), 1e-6),
      omega("omega", cells.size(), 1.0),
      mu_t("mu_t", cells.size(), 0.0),
      wallDistance("wallDistance", cells.size(), 1.0),
      wallShearStress("wallShearStress", cells.size(), 0.0),
      grad_k("grad_k", cells.size(), Vector(0.0, 0.0, 0.0)),
      grad_omega("grad_omega", cells.size(), Vector(0.0, 0.0, 0.0)),
      allFaces(faces),
      allCells(cells),
      bcManager(bc),
      gradientScheme(gradientScheme),
      F1("F1", cells.size(), 0.0),
      F2("F2", cells.size(), 0.0),
      production_k("production_k", cells.size(), 0.0),
      production_omega("production_omega", cells.size(), 0.0),
      cross_diffusion("cross_diffusion", cells.size(), 0.0),
      enableTransient(false),
      dt(0.001),
      theta(0.5),
      k_old("k_old", cells.size(), 1e-6),
      omega_old("omega_old", cells.size(), 1.0),
      rho(1.225)  // Default air density
{
    // Initialize matrix constructor for equation solving
    matrixConstruct = std::make_unique<Matrix>
    (
        allFaces,
        allCells,
        bcManager
    );
    
    std::cout   << "k-omega SST turbulence model initialized with "
                << allCells.size() << " cells." << std::endl;
}

KOmegaSST::~KOmegaSST() = default;

void KOmegaSST::initialize
(
    const VectorField& U_field,
    Scalar nu_lam
)
{
    std::cout   << "\n=== Initializing k-omega SST Turbulence Model ==="
                << std::endl;
    
    // Step 1: Calculate wall distance using Poisson equation
    calculateWallDistance();
    
    // Step 2: Initialize turbulence fields based on flow conditions
    Scalar I_turbulence = 0.05;  // 5% turbulence intensity
    Scalar l_turbulent = 0.1;    // Turbulent length scale (adjust based on geometry)
    
    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        Scalar U_mag = U_field[i].magnitude();
        
        // Initialize k based on turbulence intensity: k = 1.5 * (I * U)²
        k[i] = std::max
        (
            1.5 * (I_turbulence * U_mag) * (I_turbulence * U_mag),
            1e-8
        );
        
        // Initialize omega based on length scale: ω = k^0.5 / (C_μ^0.25 * l)
        omega[i] = std::max
        (
            std::sqrt(k[i]) / (std::pow(constants.C_mu, 0.25) * l_turbulent),
            1e-4
        );
        
        // Initialize turbulent kinematic viscosity: νₜ = k / ω  
        mu_t[i] = k[i] / omega[i];  // This is actually nu_t (kinematic)
    }
    
    // Step 3: Apply turbulence boundary conditions
    applyTurbulenceBoundaryConditions("k", k, U_field, nu_lam);
    applyTurbulenceBoundaryConditions("omega", omega, U_field, nu_lam);
    
    std::cout   << "Turbulence fields initialized successfully." << std::endl;
}

void KOmegaSST::calculateWallDistance()
{
    std::cout   << "  Computing wall distance using Poisson equation..." 
                << std::endl;
    
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
    VectorField zero_velocity("zero_velocity", allCells.size(), Vector(0.0, 0.0, 0.0));
    ScalarField phi_wall_source("phi_wall_source", allCells.size(), 0.0);
    ScalarField Gamma_phi_wall("Gamma_phi_wall", allCells.size(), 1.0);
    UpwindScheme upwds;
    
    matrixConstruct->buildMatrix
    (
        phi,
        phi_wall_source, // Source term (zero for wall distance)
        zero_velocity,   // No convection
        Gamma_phi_wall,  // Unit diffusion coefficient per cell
        upwds,
        gradientScheme,
        "phi_wall"
    );
    
    // Modify source term: b = -1 (Laplacian = -1)
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>
    (
        matrixConstruct->getMatrixA()
    );

    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>
    (
        matrixConstruct->getVectorB()
    );
    
    // Set source term to -1 for all cells
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        b_vector(i) = -allCells[i].volume();  // -1 * cell_volume
    }
    
    // Apply boundary condition: φ = 0 at walls only; Neumann elsewhere
    // Build face-to-patch map
    std::map<size_t, const BoundaryPatch*> faceToPatch;

    for (const auto& patch : bcManager.patches())
    {
        for (size_t i = patch.firstFaceIndex(); i <= patch.lastFaceIndex(); ++i)
        {
            faceToPatch[i] = &patch;
        }
    }

    for (const auto& face : allFaces)
    {
        if (!face.isBoundary()) continue;

        size_t P = face.ownerCell();

        const BoundaryPatch* patch = faceToPatch.at(face.id());

        if (patch->type() == BoundaryConditionType::WALL)
        {
            A_matrix.coeffRef(P, P) += 1e12;  // Dirichlet pin
            b_vector(P) = 0.0;               // φ = 0 at wall
        }
        else
        {
            // Zero normal gradient elsewhere: no extra constraint (already implicit)
        }
    }
    
    // Solve the Poisson equation
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> phi_solution(allCells.size());
    phi_solution.setZero();
    
    bool solved = LinearSolvers::BiCGSTAB
    (
        phi_solution,
        A_matrix,
        b_vector,
        1e-8,
        1000,
        "phi_wall"
    );
    
    if (solved)
    {
        for (size_t i = 0; i < allCells.size(); ++i)
        {
            phi[i] = phi_solution(i);
        }
    }
    
    // Compute gradients of phi field
    VectorField grad_phi("grad_phi", allCells.size(), Vector(0.0, 0.0, 0.0));
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        grad_phi[i] = gradientScheme.CellGradient(i, phi, allCells);
    }
    
    // Compute wall distance: d = √(|∇φ|² + 2φ) - |∇φ|
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        Scalar grad_phi_mag = grad_phi[i].magnitude();
        Scalar discriminant = grad_phi_mag * grad_phi_mag + 2.0 * phi[i];
        
        if (discriminant > 0)
        {
            wallDistance[i] = std::sqrt(discriminant) - grad_phi_mag;
        }
        else
        {
            wallDistance[i] = 1e-6;  // Minimum wall distance
        }
    
        // Ensure positive wall distance with a much smaller floor
        wallDistance[i] = std::max(wallDistance[i], 1e-9);
    }
    
    std::cout << "  Wall distance calculation completed." << std::endl;
}

// Main solve method that calculates gradients internally
void KOmegaSST::solve
(
    const VectorField& U_field,
    Scalar nu_lam
)
{
    // Extract velocity components
    ScalarField Ux("Ux", allCells.size());
    ScalarField Uy("Uy", allCells.size());
    ScalarField Uz("Uz", allCells.size());
    
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        Ux[i] = U_field[i].x();
        Uy[i] = U_field[i].y();
        Uz[i] = U_field[i].z();
    }
    
    // Calculate velocity gradients
    VectorField gradUx("gradUx", allCells.size(), Vector(0.0, 0.0, 0.0));
    VectorField gradUy("gradUy", allCells.size(), Vector(0.0, 0.0, 0.0));
    VectorField gradUz("gradUz", allCells.size(), Vector(0.0, 0.0, 0.0));
    
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        gradUx[i] = gradientScheme.CellGradient(i, Ux, allCells);
        gradUy[i] = gradientScheme.CellGradient(i, Uy, allCells);
        gradUz[i] = gradientScheme.CellGradient(i, Uz, allCells);
    }
    
    std::vector<VectorField> gradUVec = { gradUx, gradUy, gradUz };
    solve(U_field, gradUVec, nu_lam);
}

// Legacy vector-based interface (unchanged)
void KOmegaSST::solve
(
    const VectorField& U_field, 
    const std::vector<VectorField>& gradU,
    Scalar nu_lam
)
{
    /**
     * Complete k-omega SST solution sequence:
     * 1. Calculate blending functions F1 and F2
     * 2. Solve omega equation with near-wall treatment
     * 3. Solve k equation
     * 4. Calculate turbulent viscosity
     * 5. Apply wall corrections
     */
    
    // Step 1: Calculate gradients of turbulence quantities
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        grad_k[i] = gradientScheme.CellGradient(i, k, allCells);
        grad_omega[i] = gradientScheme.CellGradient(i, omega, allCells);
    }
    
    // Step 2: Calculate blending functions
    calculateBlendingFunctions(gradU, nu_lam);
    
    // Step 3: Calculate production terms
    calculateProductionTerms(gradU);
    
    // Step 4: Calculate cross-diffusion term
    calculateCrossDiffusion();
    
    // Step 5: Solve omega equation
    solveOmegaEquation(U_field, gradU, nu_lam);
    
    // Step 6: Apply near-wall treatment for omega
    applyNearWallTreatmentOmega(nu_lam);
    
    // Step 7: Solve k equation
    solveKEquation(U_field, gradU, nu_lam);
    
    // Step 8: Calculate turbulent viscosity
    calculateTurbulentViscosity(U_field, gradU, nu_lam);
    
    // Step 9: Apply wall corrections
    applyWallCorrections();
    
    // Step 10: Calculate wall shear stress
    calculateWallShearStress(U_field, nu_lam);
}

void KOmegaSST::solveOmegaEquation
(
    const VectorField& U_field,
    const std::vector<VectorField>&,
    Scalar nu_lam
)
{   
    /**
     * Omega transport equation (kinematic formulation):
     * ∂ω/∂t + ∇·(Uω) = γ·S² - β·ω² + ∇·[(ν + σ_ω·νₜ)∇ω] + D_ω
     * 
     * Where:
     * - γ·S² = production term (kinematic)
     * - D_ω = cross-diffusion term (only for SST model)
     * - γ, β, σ_ω are blended constants
     */
    
    std::cout << "  Solving omega transport equation..." << std::endl;
    
    // omega_old is already stored as a member variable
    
    // Construct transport matrix for omega (variable effective diffusion: ν + σ_ω ν_t)
    // Build per-cell effective diffusion Gamma = ν_lam + σ_ω * ν_t (kinematic viscosity)
    ScalarField GammaOmega("GammaOmega", allCells.size());

    for (size_t i = 0; i < allCells.size(); ++i)
    {
        Scalar sigma_omega = 
            F1[i] * constants.sigma_omega1 
          + (S(1.0) - F1[i]) * constants.sigma_omega2;

        GammaOmega[i] = nu_lam + sigma_omega * mu_t[i];  // All kinematic
    }

    // Construct transport matrix for omega with variable diffusion
    ScalarField omega_source("omega_source", allCells.size(), 0.0);
    CentralDifferenceScheme cds_omega;
    
    matrixConstruct->buildMatrix
    (
        omega,
        omega_source, // Source term (zero for omega)
        U_field,
        GammaOmega,
        cds_omega,
        gradientScheme,
        "omega"
    );
    
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>
    (
        matrixConstruct->getMatrixA()
    );

    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>
    (
        matrixConstruct->getVectorB()
    );
    
    // Add omega-specific source terms and modify diffusion
    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        Scalar cellVolume = allCells[i].volume();
        
        // Blended constants
        const Scalar blended_gamma = 
            F1[i] * constants.gamma_1 
          + (1.0 - F1[i]) * constants.gamma_2;

        const Scalar beta = 
            F1[i] * constants.beta_1 
          + (1.0 - F1[i]) * constants.beta_2;
        
        // Production term for omega: γ * S² (kinematic formulation)
        // production_omega is already S² from calculateProductionTerms
        Scalar P_omega = blended_gamma * production_omega[i];

        b_vector(i) += P_omega * cellVolume;
        
        // Destruction term: -β·ω² (kinematic formulation)
        const Scalar destruction = beta * omega[i] * omega[i];
        A_matrix.coeffRef(i, i) += destruction / omega[i] * cellVolume;  // Linearized
        
        // Cross-diffusion term (SST specific)
        b_vector(i) += cross_diffusion[i] * cellVolume;
        
        // Enhanced diffusion due to turbulent viscosity
        // This would require modifying the matrix construction for variable diffusion
        // For simplicity, we approximate this effect
    }
    
    // Solve omega equation
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> omega_solution(allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        omega_solution(i) = omega[i];
    }
    
    bool solved = LinearSolvers::BiCGSTAB
    (
        omega_solution,
        A_matrix,
        b_vector,
        1e-8,
        1000,
        "omega"
    );
    
    if (solved)
    {
        for (size_t i = 0; i < allCells.size(); ++i)
        {
            omega[i] = std::max(omega_solution(i), 1e-6);  // Ensure positive omega
        }
    }
    
    // Apply boundary conditions
    applyTurbulenceBoundaryConditions("omega", omega, U_field, nu_lam);
}

void KOmegaSST::applyNearWallTreatmentOmega(Scalar nu_lam) 
{
    /**
     * Near-wall treatment for omega:
     * For cells very close to the wall (y+ < 2), use the viscous sublayer relation:
     * ω_wall = 6μ/(ρβ₁y²) = 6ν/(β₁y²)
     */
        
    std::cout << "  Applying near-wall treatment for omega..." << std::endl;
    
    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        Scalar y = wallDistance[i];
        // Trigger near-wall if y is very small (tolerant threshold)
        if (y < 5e-6) 
        {
            omega[i] = 6.0 * nu_lam / (constants.beta_1 * y * y + 1e-20);
            omega[i] = std::min(omega[i], 1e6);
            omega[i] = std::max(omega[i], 1e-4);
        }
    }
}

void KOmegaSST::solveKEquation
(
    const VectorField& U_field,
    const std::vector<VectorField>&,
    Scalar nu_lam
)
{
    /**
     * k transport equation (kinematic formulation):
     * ∂k/∂t + ∇·(Uk) = P_k - β*·kω + ∇·[(ν + σ_k·νₜ)∇k]
     * 
     * Where:
     * - P_k = limited production term (kinematic)
     * - β* = 0.09 (destruction coefficient) 
     * - σ_k = blended diffusion coefficient
     */
    
    std::cout << "  Solving k transport equation..." << std::endl;
    
    // k_old is already stored as a member variable
    
    // Construct transport matrix for k (variable effective diffusion: μ + σ_k μ_t)
    ScalarField GammaK("GammaK", allCells.size());

    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        Scalar sigma_k = 
            F1[i] * constants.sigma_k1 
          + (S(1.0) - F1[i]) * constants.sigma_k2;
          
        GammaK[i] = nu_lam + sigma_k * mu_t[i];  // All kinematic
    }
    ScalarField k_source("k_source", allCells.size(), 0.0);
    CentralDifferenceScheme cds_k;
    
    matrixConstruct->buildMatrix
    (
        k,
        k_source, // Source term (zero for k)
        U_field,
        GammaK,
        cds_k,
        gradientScheme,
        "k"
    );
    
    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>
    (
        matrixConstruct->getMatrixA()
    );

    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>
    (
        matrixConstruct->getVectorB()
    );
    
    // Add k-specific source terms
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        Scalar cellVolume = allCells[i].volume();
        
        // Production term (limited to prevent unrealistic values)
        const Scalar P_k_limited = limitProduction(production_k[i], i);
        b_vector(i) += P_k_limited * cellVolume;
        
        // Destruction term: -β*·kω (kinematic formulation)  
        const Scalar destruction = constants.beta_star * omega[i];
        A_matrix.coeffRef(i, i) += destruction * cellVolume;
    }
    
    // Solve k equation
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> k_solution(allCells.size());
    for (size_t i = 0; i < allCells.size(); ++i) {
        k_solution(i) = k[i];
    }
    
    bool solved = LinearSolvers::BiCGSTAB
    (
        k_solution,
        A_matrix,
        b_vector,
        1e-8,
        1000,
        "k"
    );
    
    if (solved) 
    {
        for (size_t i = 0; i < allCells.size(); ++i) 
        {
            k[i] = std::max(k_solution(i), 1e-8);  // Ensure positive k
        }
    }
    
    // Apply boundary conditions
    applyTurbulenceBoundaryConditions("k", k, U_field, nu_lam);
}

void KOmegaSST::calculateTurbulentViscosity
(
    const VectorField&,
    const std::vector<VectorField>& gradU,
    Scalar nu_lam
)
{
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
    
    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        // Calculate strain rate magnitude S (not vorticity!)
        Scalar S = 0.0;
        if (gradU.size() >= 3) 
        {  // Ensure we have all velocity gradients
            // S = √(2SᵢⱼSᵢⱼ) where Sᵢⱼ is the strain rate tensor
            Scalar S11 = gradU[0][i].x();                                    // ∂u/∂x
            Scalar S22 = gradU[1][i].y();                                    // ∂v/∂y  
            Scalar S33 = gradU[2][i].z();                                    // ∂w/∂z
            Scalar S12 = 0.5 * (gradU[0][i].y() + gradU[1][i].x());          // 0.5(∂u/∂y + ∂v/∂x)
            Scalar S13 = 0.5 * (gradU[0][i].z() + gradU[2][i].x());          // 0.5(∂u/∂z + ∂w/∂x)
            Scalar S23 = 0.5 * (gradU[1][i].z() + gradU[2][i].y());          // 0.5(∂v/∂z + ∂w/∂y)
            
            S = std::sqrt(
                2.0 * (S11*S11 + S22*S22 + S33*S33 
              + 2.0*(S12*S12 + S13*S13 + S23*S23))
            );
        }
        else
        {
            // Simplified 2D case
            Scalar S11 = gradU[0][i].x();
            Scalar S22 = gradU[1][i].y();
            Scalar S12 = 0.5 * (gradU[0][i].y() + gradU[1][i].x());
            
            S = std::sqrt(2.0 * (S11*S11 + S22*S22 + 2.0*S12*S12));
        }
        
        // SST turbulent viscosity formula: νₜ = a₁ * k / max(a₁ * ω, S * F₂)
        Scalar denominator = std::max(constants.a1 * omega[i], S * F2[i]);
        mu_t[i] = constants.a1 * k[i] / (denominator + 1e-20);  // This is actually nu_t (kinematic)
        
        // Limit turbulent viscosity to reasonable values
        mu_t[i] = std::min(mu_t[i], 1000.0 * nu_lam);  // Max 1000 times laminar kinematic viscosity
        mu_t[i] = std::max(mu_t[i], 0.0);              // Ensure non-negative
    }
}

void KOmegaSST::applyWallCorrections() 
{
    /**
     * Wall corrections for turbulent viscosity:
     * - μₜ = 0 at walls (no-slip condition)
     * - Smooth transition from wall to freestream
     * - Apply damping functions based on y+ or wall distance
     */
    
    std::cout << "  Applying wall corrections..." << std::endl;
    
    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        Scalar y = wallDistance[i];
        
        // Wall damping function (exponential decay near wall)
        // f_μ = 1 - exp(-y/δ) where δ is characteristic length
        Scalar delta = 1e-4;  // Characteristic damping length
        Scalar f_mu = 1.0 - std::exp(-y / delta);
        
        // Apply damping to turbulent viscosity
        mu_t[i] *= f_mu;
        
        // Ensure μₜ = 0 at walls (y → 0)
        if (y < 1e-6) 
        {
            mu_t[i] = 0.0;
        }
    }
}

void KOmegaSST::calculateWallShearStress(const VectorField& U_field, Scalar nu_lam)
{
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
    
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        Scalar y = wallDistance[i];
        
        if (y < 1e-3)   // Near-wall cells
        {
            // Approximate wall shear using first cell velocity and distance
            Scalar U_magnitude = U_field[i].magnitude();
            
            // Wall shear stress: τ = μ * ∂U/∂y ≈ μ * U / y = ρ * ν * U / y
            wallShearStress[i] = rho * nu_lam * U_magnitude / (y + 1e-12);
            
            // Limit to reasonable values
            wallShearStress[i] = std::min(wallShearStress[i], 1000.0);
        } 
        else 
        {
            wallShearStress[i] = 0.0;  // Far from wall
        }
    }
}

ScalarField KOmegaSST::getEffectiveViscosity(Scalar nu_lam) const
{
    /**
     * Calculate effective kinematic viscosity: ν_eff = ν_lam + νₜ
     * This is used in momentum equations for turbulent flow
     */
    
    ScalarField nu_eff("nu_eff", allCells.size());
    
    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        nu_eff[i] = nu_lam + mu_t[i];  // mu_t is actually nu_t (kinematic)
    }
    
    return nu_eff;
}

// Private helper methods implementation

void KOmegaSST::calculateBlendingFunctions
(
    const std::vector<VectorField>&,
    Scalar nu_lam
)
{
    /**
     * Calculate SST blending functions F1 and F2:
     * 
     * F1 = tanh(arg1⁴) where arg1 = min(max(√k/(β*ωy), 500μ/(ρωy²)), 4ρσω2k/(CDkωy²))
     * F2 = tanh(arg2²) where arg2 = max(2√k/(β*ωy), 500μ/(ρωy²))
     * 
     * These functions blend between k-ω (F1=1) and k-ε (F1=0) formulations
     */
    
    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        Scalar y = wallDistance[i];
        Scalar sqrt_k = std::sqrt(k[i]);
        
        // Avoid division by zero
        y = std::max(y, 1e-12);
        omega[i] = std::max(omega[i], 1e-12);
        
        // Arguments for blending functions
        Scalar arg1_1 = sqrt_k / (constants.beta_star * omega[i] * y);
        Scalar arg1_2 = 500.0 * rho * nu_lam / (omega[i] * y * y);
        
        // Calculate CDkw for arg1_3 (cross-diffusion for blending function)
        Scalar dot_product = dot(grad_k[i], grad_omega[i]);

        Scalar CDkw = 
            std::max
            (
                2.0 * constants.sigma_omega2 / omega[i] * dot_product,
                1e-10
            );

        Scalar arg1_3 = 
            4.0 * constants.sigma_omega2 * k[i] / (CDkw * y * y);
        
        Scalar arg1 = std::min(std::max(arg1_1, arg1_2), arg1_3);
        F1[i] = std::tanh(std::pow(arg1, 4.0));
        
        Scalar arg2 = 
            std::max
            (
                2.0 * sqrt_k / (constants.beta_star * omega[i] * y),
                arg1_2
            );

        F2[i] = std::tanh(arg2 * arg2);
    }
}

void KOmegaSST::calculateProductionTerms
(
    const std::vector<VectorField>& gradU
)
{
    /**
     * Calculate production terms for k and omega:
     * 
     * P_k = μₜ * S² where S is strain rate magnitude
     * P_ω = (γ/νₜ) * P_k where γ is blended coefficient
     */
    
    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        // Calculate strain rate magnitude S = √(2SᵢⱼSᵢⱼ)
        Scalar S = 0.0;
        if (gradU.size() >= 3)
        {
            // Full 3D strain rate tensor
            Scalar S11 = gradU[0][i].x();                                    // ∂u/∂x
            Scalar S22 = gradU[1][i].y();                                    // ∂v/∂y  
            Scalar S33 = gradU[2][i].z();                                    // ∂w/∂z
            Scalar S12 = 0.5 * (gradU[0][i].y() + gradU[1][i].x());          // 0.5(∂u/∂y + ∂v/∂x)
            Scalar S13 = 0.5 * (gradU[0][i].z() + gradU[2][i].x());          // 0.5(∂u/∂z + ∂w/∂x)
            Scalar S23 = 0.5 * (gradU[1][i].z() + gradU[2][i].y());          // 0.5(∂v/∂z + ∂w/∂y)
            
            S = std::sqrt
                (
                    2.0 * (S11*S11 + S22*S22 + S33*S33 
                  + 2.0*(S12*S12 + S13*S13 + S23*S23))
                );
        }
        else
        {
            // Simplified 2D case
            Scalar S11 = gradU[0][i].x();
            Scalar S22 = gradU[1][i].y();
            Scalar S12 = 0.5 * (gradU[0][i].y() + gradU[1][i].x());
            
            S = std::sqrt(2.0 * (S11*S11 + S22*S22 + 2.0*S12*S12));
        }
        
        // Production of k: P_k = νₜ * S² (kinematic formulation)
        production_k[i] = mu_t[i] * S * S;  // mu_t is actually nu_t
        
        // Production of omega: P_ω = (γ/νₜ) * P_k = γ * S²
        Scalar gamma = 
            F1[i] * constants.gamma_1 + (1.0 - F1[i]) * constants.gamma_2;

        production_omega[i] = gamma * S * S;
    }
}

void KOmegaSST::calculateCrossDiffusion()
{
    /**
     * Calculate cross-diffusion term for SST model (kinematic formulation):
     * D_ω = 2(1-F1) * σ_ω2 / ω * ∇k · ∇ω
     * 
     * This term appears only in the omega equation and is crucial for SST behavior
     */
        
    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        // Cross-diffusion: 2(1-F1) * σ_ω2 / ω * ∇k · ∇ω (kinematic formulation)
        Scalar dot_product = dot(grad_k[i], grad_omega[i]);
        
        cross_diffusion[i] = 
            2.0 * (1.0 - F1[i]) * constants.sigma_omega2
          / (omega[i] + 1e-12) * dot_product;
        
        // Ensure positive contribution (only for SST formulation)
        cross_diffusion[i] = std::max(cross_diffusion[i], 0.0);
    }
}

void KOmegaSST::applyTurbulenceBoundaryConditions
(
    const std::string& fieldName,
    ScalarField& field,
    const VectorField&,
    Scalar nu_lam
)
{
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
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        if (wallDistance[i] < 1e-6)   // At wall
        {
            if (fieldName == "k")
            {
                field[i] = 0.0;  // k = 0 at wall
            } 
            else if (fieldName == "omega") 
            {
                // ω_wall = 6μ/(ρβ₁y²) for viscous sublayer
                Scalar y = std::max(wallDistance[i], 1e-8);
                field[i] = 6.0 * nu_lam / (constants.beta_1 * y * y);
                field[i] = std::min(field[i], 1e6);  // Limit maximum value
            }
        }
    }
}

Scalar KOmegaSST::calculateYPlus
(
    size_t cellIdx,
    const VectorField&,
    Scalar nu_lam
) const
{
    /**
     * Calculate y+ = ρ * u_τ * y / μ
     * where u_τ = √(τ_wall/ρ) is friction velocity
     */
    
    Scalar y = wallDistance[cellIdx];
    Scalar tau_wall = wallShearStress[cellIdx];
    Scalar u_tau = std::sqrt(tau_wall);
    
    return u_tau * y / nu_lam;
}

Scalar KOmegaSST::limitProduction(Scalar P_k, size_t cellIdx) const
{
    /**
     * Limit production to prevent unrealistic values:
     * P_k = min(P_k, 10 * β* * k * ω)  (kinematic formulation)
     * 
     * This prevents excessive production in stagnation regions
     */
    
    Scalar limit = 
        10.0 * constants.beta_star * k[cellIdx] * omega[cellIdx];

    return std::min(P_k, limit);
}

void KOmegaSST::setModelConstants(Scalar C_mu, Scalar beta_1, Scalar beta_2)
{
    constants.C_mu = C_mu;
    constants.beta_1 = beta_1;
    constants.beta_2 = beta_2;
}

void KOmegaSST::setTransientMode(bool enable, Scalar dt_new, Scalar theta_new)
{
    enableTransient = enable;
    dt = dt_new;
    theta = theta_new;
}

void KOmegaSST::updateOldValues() 
{
    // Store current values as old values for next time step
    k_old = k;
    omega_old = omega;
} 