/******************************************************************************
 * @file KOmegaSST.cpp
 * @brief Implementation of k-omega SST turbulence model
 *****************************************************************************/

#include <cmath>
#include <algorithm>
#include <iostream>
#include <set>

#include "kOmegaSST.hpp"
#include "Matrix.hpp"
#include "ConvectionScheme.hpp"


// ************************ Constructor & Destructor ************************

KOmegaSST::KOmegaSST
(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const BoundaryConditions& bc,
    const GradientScheme& gradientScheme,
    const ConvectionScheme& kScheme,
    const ConvectionScheme& omegaScheme
)
    : k_("k", cells.size(), 1e-6),
      omega_("omega", cells.size(), 1.0),
      nu_t_("nu_t", cells.size(), 0.0),
      wallDistance_("wallDistance", cells.size(), 1.0),
      wallShearStress_("wallShearStress", cells.size(), 0.0),
      grad_k_("grad_k", cells.size(), Vector(0.0, 0.0, 0.0)),
      grad_omega_("grad_omega", cells.size(), Vector(0.0, 0.0, 0.0)),
      allFaces_(faces),
      allCells_(cells),
      bcManager_(bc),
      gradientScheme_(gradientScheme),
      kConvectionScheme_(kScheme),
      omegaConvectionScheme_(omegaScheme),
      F1_("F1", cells.size(), 0.0),
      F2_("F2", cells.size(), 0.0),
      production_k_("production_k", cells.size(), 0.0),
      production_omega_("production_omega", cells.size(), 0.0),
      cross_diffusion_("cross_diffusion", cells.size(), 0.0),
      rho_(1.225)  // Default air density
{
    // Initialize matrix constructor for equation solving
    matrixConstruct_ = std::make_unique<Matrix>
    (
        allFaces_,
        allCells_,
        bcManager_
    );
    
    std::cout   << "k-omega SST turbulence model initialized with "
                << allCells_.size() << " cells." << std::endl;
}

KOmegaSST::~KOmegaSST() = default;


// ******************************* Public Methods ******************************

void KOmegaSST::initialize
(
    const VectorField& U_field,
    Scalar nu_lam
)
{
    std::cout   << "\n=== Initializing k-omega SST Turbulence Model ==="
                << std::endl;
    
    calculateWallDistance();
    
    Scalar I_turbulence = 0.05;  // 5% turbulence intensity
    Scalar l_turbulent = 0.1;    // Turbulent length scale
    
    for (size_t i = 0; i < allCells_.size(); ++i) 
    {
        Scalar U_mag = U_field[i].magnitude();
        
        // Initialize k based on turbulence intensity: k = 1.5 * (I * U)²
        k_[i] = std::max
        (
            1.5 * (I_turbulence * U_mag) * (I_turbulence * U_mag),
            1e-8
        );
        
        // Initialize omega based on length scale: ω = k^0.5 / (C_μ^0.25 * l)
        omega_[i] = std::max
        (
            std::sqrt(k_[i]) / (std::pow(constants_.C_mu, 0.25) * l_turbulent),
            1e-4
        );
        
        // Initialize turbulent kinematic viscosity: νₜ = k / ω
        nu_t_[i] = k_[i] / omega_[i];
    }
    
    applyTurbulenceBoundaryConditions("k", k_, U_field, nu_lam);
    applyTurbulenceBoundaryConditions("omega", omega_, U_field, nu_lam);
    
    std::cout   << "Turbulence fields initialized successfully." << std::endl;
}

void KOmegaSST::calculateWallDistance()
{
    std::cout   << "  Computing wall distance using meshWave method..."
                << std::endl;

    /**
     * meshWave (iterative propagation) approach for wall distance:
     *
     * This method propagates exact geometric distances from walls through
     * the mesh using face connectivity. It is equivalent to OpenFOAM's
     * meshWave method and provides exact distances on regular meshes.
     *
     * Algorithm:
     * 1. Initialize wall-adjacent cells with exact distance to wall face
     * 2. Iteratively propagate distances through internal faces
     * 3. Each cell's distance = min(current, neighbor_dist + cell_spacing)
     * 4. Converges when no distances change (typically 30-50 iterations)
     */

    size_t numCells = allCells_.size();

    // Step 1: Initialize all cell distances to large value
    wallDistance_.setAll(1e10);

    // Step 2: Initialize wall-adjacent cells with exact geometric distance
    const auto& patchMap = bcManager_.faceToPatchMap();

    for (const auto& face : allFaces_)
    {
        if (!face.isBoundary()) continue;

        // O(1) patch lookup using cached map
        auto it = patchMap.find(face.idx());
        const BoundaryPatch* patch = (it != patchMap.end()) ? it->second : nullptr;

        if (!patch || patch->type() != BoundaryConditionType::WALL) continue;

        size_t cellIdx = face.ownerCell();
        Vector cellCenter = allCells_[cellIdx].centroid();
        Vector faceCenter = face.centroid();
        Vector normal = face.normal();  // Already unit vector

        // Exact perpendicular distance from cell center to wall face
        Vector cellToFace = faceCenter - cellCenter;
        Scalar dist = std::abs(dot(cellToFace, normal));

        // Keep minimum distance if cell touches multiple walls
        wallDistance_[cellIdx] = std::min(wallDistance_[cellIdx], dist);
    }

    // Step 3: Iterative propagation through mesh
    const size_t maxIterations = 100;
    const Scalar tolerance = 1e-12;

    for (size_t iter = 0; iter < maxIterations; ++iter)
    {
        Scalar maxChange = 0.0;

        // Iterate through all internal faces
        for (const auto& face : allFaces_)
        {
            if (face.isBoundary()) continue;

            size_t owner = face.ownerCell();
            size_t neighbor = face.neighborCell().value();

            Vector ownerCenter = allCells_[owner].centroid();
            Vector neighborCenter = allCells_[neighbor].centroid();
            Scalar cellDist = (ownerCenter - neighborCenter).magnitude();

            // Update owner from neighbor
            Scalar newDistOwner = wallDistance_[neighbor] + cellDist;
            if (newDistOwner < wallDistance_[owner])
            {
                Scalar change = wallDistance_[owner] - newDistOwner;
                maxChange = std::max(maxChange, change);
                wallDistance_[owner] = newDistOwner;
            }

            // Update neighbor from owner
            Scalar newDistNeighbor = wallDistance_[owner] + cellDist;
            if (newDistNeighbor < wallDistance_[neighbor])
            {
                Scalar change = wallDistance_[neighbor] - newDistNeighbor;
                maxChange = std::max(maxChange, change);
                wallDistance_[neighbor] = newDistNeighbor;
            }
        }

        // Check convergence
        if (maxChange < tolerance)
        {
            break;
        }

        if (iter == maxIterations - 1)
        {
            std::cout   << "  WARNING: Reached max iterations " 
                        << "without convergence"
                        << std::endl;
        }
    }

    // Step 4: Apply minimum threshold and compute statistics
    Scalar minDist = 1e10;
    Scalar maxDist = 0.0;

    for (size_t i = 0; i < numCells; ++i)
    {
        // Apply minimum threshold to prevent numerical issues
        wallDistance_[i] = std::max(wallDistance_[i], 1e-12);

        minDist = std::min(minDist, wallDistance_[i]);
        maxDist = std::max(maxDist, wallDistance_[i]);
    }

    std::cout << "  Wall distance calculation completed." << std::endl;
}

// Main solve method that calculates gradients internally
void KOmegaSST::solve
(
    const VectorField& U_field,
    const FaceFluxField& flowRateFace,
    Scalar nu_lam,
    Scalar alpha_k,
    Scalar alpha_omega
)
{
    // Extract velocity components
    ScalarField Ux("U_x", allCells_.size());
    ScalarField Uy("U_y", allCells_.size());
    ScalarField Uz("U_z", allCells_.size());

    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        Ux[i] = U_field[i].x();
        Uy[i] = U_field[i].y();
        Uz[i] = U_field[i].z();
    }

    // Calculate velocity gradients
    VectorField gradUx("gradUx", allCells_.size(), Vector(0.0, 0.0, 0.0));
    VectorField gradUy("gradUy", allCells_.size(), Vector(0.0, 0.0, 0.0));
    VectorField gradUz("gradUz", allCells_.size(), Vector(0.0, 0.0, 0.0));

    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        gradUx[i] = gradientScheme_.CellGradient(
            i, Ux, allCells_, allFaces_, bcManager_, "U_x");
        gradUy[i] = gradientScheme_.CellGradient(
            i, Uy, allCells_, allFaces_, bcManager_, "U_y");
        gradUz[i] = gradientScheme_.CellGradient(
            i, Uz, allCells_, allFaces_, bcManager_, "U_z");
    }

    // Package gradients for use in other methods
    std::vector<VectorField> gradU = {gradUx, gradUy, gradUz};

    // Calculate gradients of turbulence fields
    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        grad_k_[i] = gradientScheme_.CellGradient(
            i, k_, allCells_, allFaces_, bcManager_, "k");
        grad_omega_[i] = gradientScheme_.CellGradient(
            i, omega_, allCells_, allFaces_, bcManager_, "omega");
    }

    // Calculate blending functions F1 and F2
    calculateBlendingFunctions(gradU, nu_lam);

    // Calculate production terms
    calculateProductionTerms(gradU);

    // Calculate cross-diffusion term for SST
    calculateCrossDiffusion();

    // Solve k equation with under-relaxation
    solveKEquation(U_field, flowRateFace, gradU, nu_lam, alpha_k);

    // Solve omega equation with under-relaxation
    solveOmegaEquation(U_field, flowRateFace, gradU, nu_lam, alpha_omega);

    // Apply near-wall treatment for omega
    applyNearWallTreatmentOmega(nu_lam);

    // Calculate turbulent viscosity
    calculateTurbulentViscosity(U_field, gradU, nu_lam);

    // Apply wall corrections
    applyWallCorrections();

    // Calculate wall shear stress
    calculateWallShearStress(U_field, nu_lam);
}

void KOmegaSST::solveOmegaEquation
(
    const VectorField& U_field,
    const FaceFluxField& flowRateFace,
    const std::vector<VectorField>&,
    Scalar nu_lam,
    Scalar alpha_omega
)
{
    /**
     * Omega transport equation (kinematic formulation):
     * ∂ω/∂t + ∇·(Uω) = γ·S² - β·ω² + ∇·[(ν + σ_ω·νₜ)∇ω] + D_ω
     *
     * Where:
     * - γ·S² = production term (kinematic)
     * - D_ω = cross-diffusion term (only for SST model)
     * - γ, β, σ_ω are blended constants_
     */

    std::cout   << "  Solving omega transport equation (relaxation = "
                << alpha_omega << ")..." << std::endl;


    // Construct transport matrix for omega 
    // (variable effective diffusion: ν + σ_ω ν_t)
    // Build per-cell effective diffusion Gamma = ν_lam + σ_ω * ν_t
    ScalarField GammaOmega("GammaOmega", allCells_.size());

    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        Scalar sigma_omega =
            F1_[i] * constants_.sigma_omega1
          + (S(1.0) - F1_[i]) * constants_.sigma_omega2;

        GammaOmega[i] = nu_lam + sigma_omega * nu_t_[i];
    }

    // Construct transport matrix for omega with variable diffusion
    ScalarField omega_source("omega_source", allCells_.size(), 0.0);

    matrixConstruct_->buildMatrix
    (
        omega_,
        omega_source, // Source term (zero for omega)
        flowRateFace,
        GammaOmega,
        omegaConvectionScheme_,
        gradientScheme_,
        "omega"
    );

    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>
    (
        matrixConstruct_->getMatrixA()
    );

    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>
    (
        matrixConstruct_->getVectorB()
    );

    // Add omega-specific source terms and modify diffusion
    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        Scalar cellVolume = allCells_[i].volume();

        // Blended constants_
        const Scalar blended_gamma =
            F1_[i] * constants_.gamma_1
          + (1.0 - F1_[i]) * constants_.gamma_2;

        const Scalar beta =
            F1_[i] * constants_.beta_1
          + (1.0 - F1_[i]) * constants_.beta_2;

        // Production term for omega: γ * S² (kinematic formulation)
        // production_omega is already S² from calculateProductionTerms
        Scalar P_omega = blended_gamma * production_omega_[i];

        b_vector(i) += P_omega * cellVolume;

        // Destruction term: -β·ω² (kinematic formulation)
        const Scalar destruction = beta * omega_[i] * omega_[i];
        A_matrix.coeffRef(i, i) += destruction / omega_[i] * cellVolume;

        // Cross-diffusion term (SST specific)
        b_vector(i) += cross_diffusion_[i] * cellVolume;
    }

    // Apply implicit under-relaxation (Patankar's method)
    // Modify matrix: A' = A/α, b' = b + (1-α)/α * A_diag * φ_prev
    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        Scalar A_diag = A_matrix.coeff(i, i);
        A_matrix.coeffRef(i, i) = A_diag / alpha_omega;

        b_vector(i) +=
            (S(1.0) - alpha_omega) / alpha_omega * A_diag * omega_[i];
    }

    // Solve omega equation
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> omega_solution(allCells_.size());
    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        omega_solution(i) = omega_[i];
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
        for (size_t i = 0; i < allCells_.size(); ++i)
        {
            omega_[i] = std::max(omega_solution(i), 1e-6);  // Ensure +ve omega
        }
    }

    // Apply boundary conditions
    applyTurbulenceBoundaryConditions("omega", omega_, U_field, nu_lam);
}

void KOmegaSST::applyNearWallTreatmentOmega(Scalar nu_lam) 
{
    /**
     * Near-wall treatment for omega:
     * For (y+ < 2), use the viscous sublayer relation:
     * ω_wall = 6μ/(ρβ₁y²) = 6ν/(β₁y²)
     */
        
    std::cout << "  Applying near-wall treatment for omega..." << std::endl;
    
    for (size_t i = 0; i < allCells_.size(); ++i) 
    {
        Scalar y = wallDistance_[i];
        // Trigger near-wall if y is very small (tolerant threshold)
        if (y < 5e-6) 
        {
            omega_[i] = 6.0 * nu_lam / (constants_.beta_1 * y * y + 1e-20);
            omega_[i] = std::min(omega_[i], 1e6);
            omega_[i] = std::max(omega_[i], 1e-4);
        }
    }
}

void KOmegaSST::solveKEquation
(
    const VectorField& U_field,
    const FaceFluxField& flowRateFace,
    const std::vector<VectorField>&,
    Scalar nu_lam,
    Scalar alpha_k
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

    std::cout   << "  Solving k transport equation (relaxation = " 
                << alpha_k << ")..." << std::endl;


    // Construct transport matrix for k 
    // (variable effective diffusion: μ + σ_k μ_t)
    ScalarField GammaK("GammaK", allCells_.size());

    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        Scalar sigma_k =
            F1_[i] * constants_.sigma_k1
          + (S(1.0) - F1_[i]) * constants_.sigma_k2;

        GammaK[i] = nu_lam + sigma_k * nu_t_[i];
    }
    ScalarField k_source("k_source", allCells_.size(), 0.0);

    matrixConstruct_->buildMatrix
    (
        k_,
        k_source, // Source term (zero for k)
        flowRateFace,
        GammaK,
        kConvectionScheme_,
        gradientScheme_,
        "k"
    );

    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>
    (
        matrixConstruct_->getMatrixA()
    );

    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>
    (
        matrixConstruct_->getVectorB()
    );

    // Add k-specific source terms
    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        Scalar cellVolume = allCells_[i].volume();

        // Production term (limited to prevent unrealistic values)
        const Scalar P_k_limited = limitProduction(production_k_[i], i);
        b_vector(i) += P_k_limited * cellVolume;

        // Destruction term: -β*·kω (kinematic formulation)
        const Scalar destruction = constants_.beta_star * omega_[i];
        A_matrix.coeffRef(i, i) += destruction * cellVolume;
    }

    // Apply implicit under-relaxation (Patankar's method)
    // Modify matrix: A' = A/α, b' = b + (1-α)/α * A_diag * φ_prev
    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        Scalar A_diag = A_matrix.coeff(i, i);
        A_matrix.coeffRef(i, i) = A_diag / alpha_k;
        b_vector(i) += (S(1.0) - alpha_k) / alpha_k * A_diag * k_[i];
    }

    // Solve k equation
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> k_solution(allCells_.size());
    for (size_t i = 0; i < allCells_.size(); ++i) {
        k_solution(i) = k_[i];
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
        for (size_t i = 0; i < allCells_.size(); ++i)
        {
            k_[i] = std::max(k_solution(i), 1e-8);  // Ensure positive k
        }
    }

    // Apply boundary conditions
    applyTurbulenceBoundaryConditions("k", k_, U_field, nu_lam);
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
    
    for (size_t i = 0; i < allCells_.size(); ++i) 
    {
        // Calculate strain rate magnitude S
        Scalar S = 0.0;
        if (gradU.size() >= 3) 
        {  // Ensure we have all velocity gradients
            // S = √(2SᵢⱼSᵢⱼ) where Sᵢⱼ is the strain rate tensor
            Scalar S_11 = gradU[0][i].x();
            Scalar S_22 = gradU[1][i].y();
            Scalar S_33 = gradU[2][i].z();
            Scalar S_12 = 0.5 * (gradU[0][i].y() + gradU[1][i].x());
            Scalar S_13 = 0.5 * (gradU[0][i].z() + gradU[2][i].x());
            Scalar S_23 = 0.5 * (gradU[1][i].z() + gradU[2][i].y());
            
            S = std::sqrt
                (
                    2.0 * (S_11 * S_11 + S_22 * S_22 + S_33 * S_33 
                  + 2.0 * (S_12 * S_12 + S_13 * S_13 + S_23 * S_23))
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
        Scalar denominator = std::max(constants_.a1 * omega_[i], S * F2_[i]);
        nu_t_[i] = constants_.a1 * k_[i] / (denominator + 1e-20);
        
        // Limit turbulent viscosity to reasonable values
        nu_t_[i] = std::min(nu_t_[i], 1000.0 * nu_lam);  // Max 1000 νₜ/ν
        nu_t_[i] = std::max(nu_t_[i], 0.0);              // Ensure non-negative
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
    
    for (size_t i = 0; i < allCells_.size(); ++i) 
    {
        Scalar y = wallDistance_[i];
        
        // Wall damping function (exponential decay near wall)
        // f_μ = 1 - exp(-y/δ) where δ is characteristic length
        Scalar delta = 1e-4;  // Characteristic damping length
        Scalar f_mu = 1.0 - std::exp(-y / delta);
        
        // Apply damping to turbulent viscosity
        nu_t_[i] *= f_mu;
        
        // Ensure μₜ = 0 at walls (y → 0)
        if (y < 1e-6) 
        {
            nu_t_[i] = 0.0;
        }
    }
}

void KOmegaSST::calculateWallShearStress
(
    const VectorField& U_field,
    Scalar nu_lam
)
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
    
    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        Scalar y = wallDistance_[i];
        
        if (y < 1e-3)   // Near-wall cells
        {
            // Approximate wall shear using first cell velocity and distance
            Scalar U_magnitude = U_field[i].magnitude();
            
            // Wall shear stress: τ = μ * ∂U/∂y ≈ μ * U / y = ρ * ν * U / y
            wallShearStress_[i] = rho_ * nu_lam * U_magnitude / (y + 1e-12);
            
            // Limit to reasonable values
            wallShearStress_[i] = std::min(wallShearStress_[i], 1000.0);
        } 
        else 
        {
            wallShearStress_[i] = 0.0;  // Far from wall
        }
    }
}

ScalarField KOmegaSST::getEffectiveViscosity(Scalar nu_lam) const
{
    /**
     * Calculate effective kinematic viscosity: ν_eff = ν_lam + νₜ
     * This is used in momentum equations for turbulent flow
     */
    
    ScalarField nu_eff("nu_eff", allCells_.size());
    
    for (size_t i = 0; i < allCells_.size(); ++i) 
    {
        nu_eff[i] = nu_lam + nu_t_[i];
    }
    
    return nu_eff;
}


// ****************************** Private Methods ******************************

void KOmegaSST::calculateBlendingFunctions
(
    const std::vector<VectorField>&,
    Scalar nu_lam
)
{
    /**
     * Calculate SST blending functions F1 and F2:
     * 
     * F1 = tanh(arg1⁴) 
     * where arg1 = min(max(√k/(β*ωy), 500μ/(ρωy²)), 4ρσω2k/(CDkωy²))
     * 
     * F2 = tanh(arg2²) 
     * where arg2 = max(2√k/(β*ωy), 500μ/(ρωy²))
     * 
     * These functions blend between k-ω (F1=1) and k-ε (F1=0) formulations
     */
    
    for (size_t i = 0; i < allCells_.size(); ++i) 
    {
        Scalar y = wallDistance_[i];
        Scalar sqrt_k = std::sqrt(k_[i]);
        
        // Avoid division by zero
        y = std::max(y, 1e-12);
        omega_[i] = std::max(omega_[i], 1e-12);
        
        // Arguments for blending functions
        Scalar arg1_1 = sqrt_k / (constants_.beta_star * omega_[i] * y);
        Scalar arg1_2 = 500.0 * nu_lam / (omega_[i] * y * y);

        // Calculate CDkw for arg1_3 (cross-diffusion for blending function)
        Scalar dot_product = dot(grad_k_[i], grad_omega_[i]);

        Scalar CDkw =
            std::max
            (
                2.0 * constants_.sigma_omega2 / omega_[i] * dot_product,
                1e-10
            );

        Scalar arg1_3 =
            4.0 * constants_.sigma_omega2 * k_[i] / (CDkw * y * y);
        
        Scalar arg1 = std::min(std::max(arg1_1, arg1_2), arg1_3);
        F1_[i] = std::tanh(std::pow(arg1, 4.0));
        
        Scalar arg2 = 
            std::max
            (
                2.0 * sqrt_k / (constants_.beta_star * omega_[i] * y),
                arg1_2
            );

        F2_[i] = std::tanh(arg2 * arg2);
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
    
    for (size_t i = 0; i < allCells_.size(); ++i) 
    {
        // Calculate strain rate magnitude S = √(2SᵢⱼSᵢⱼ)
        Scalar S = 0.0;
        if (gradU.size() >= 3)
        {
            // Full 3D strain rate tensor
            Scalar S11 = gradU[0][i].x();
            Scalar S22 = gradU[1][i].y();
            Scalar S33 = gradU[2][i].z();
            Scalar S12 = 0.5 * (gradU[0][i].y() + gradU[1][i].x());
            Scalar S13 = 0.5 * (gradU[0][i].z() + gradU[2][i].x());
            Scalar S23 = 0.5 * (gradU[1][i].z() + gradU[2][i].y());
            
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
        production_k_[i] = nu_t_[i] * S * S;
        
        // Production of omega: P_ω = (γ/νₜ) * P_k = γ * S²
        Scalar gamma = 
            F1_[i] * constants_.gamma_1 + (1.0 - F1_[i]) * constants_.gamma_2;

        production_omega_[i] = gamma * S * S;
    }
}

void KOmegaSST::calculateCrossDiffusion()
{
    /**
     * Calculate cross-diffusion term for SST model (kinematic formulation):
     * D_ω = 2(1-F1) * σ_ω2 / ω * ∇k · ∇ω
     * 
     */
        
    for (size_t i = 0; i < allCells_.size(); ++i) 
    {
        // Cross-diffusion: 2(1-F1) * σ_ω2 / ω * ∇k · ∇ω 
        //(kinematic formulation)
        Scalar dot_product = dot(grad_k_[i], grad_omega_[i]);
        
        cross_diffusion_[i] = 
            2.0 * (1.0 - F1_[i]) * constants_.sigma_omega2
          / (omega_[i] + 1e-12) * dot_product;
        
        // Ensure positive contribution (only for SST formulation)
        cross_diffusion_[i] = std::max(cross_diffusion_[i], 0.0);
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
     * This method iterates through actual boundary faces and applies
     * the boundary conditions specified in the setup file.
     *
     * Supported BC types:
     * - fixedValue: Set owner cell value to specified BC value
     * - zeroGradient: Owner cell value unchanged (gradient = 0)
     * - Wall treatment: Special handling for k and omega at walls
     */

    // Build face-to-patch map for efficient BC lookup
    std::map<size_t, const BoundaryPatch*> faceToPatch;
    for (const auto& patch : bcManager_.patches())
    {
        for 
        (
            size_t faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx(); 
            ++faceIdx
        )
        {
            faceToPatch[faceIdx] = &patch;
        }
    }

    // Iterate through all boundary faces and apply BCs
    for (const auto& face : allFaces_)
    {
        if (!face.isBoundary()) continue;

        size_t ownerCell = face.ownerCell();
        const BoundaryPatch* patch = faceToPatch.at(face.idx());

        // Query boundary condition from bcManager_
        const BoundaryData* bc =
            bcManager_.fieldBC(patch->patchName(), fieldName);

        if (bc == nullptr)
        {
            // No BC specified for this field on this patch
            continue;
        }

        // Apply boundary condition based on type
        switch (bc->type())
        {
            case BCType::FIXED_VALUE:
            {
                // Set owner cell to fixed boundary value
                field[ownerCell] = bc->scalarValue();
                break;
            }

            case BCType::ZERO_GRADIENT:
            {
                // Zero gradient: owner cell value remains unchanged
                // (interior solution determines the value)
                break;
            }

            case BCType::FIXED_GRADIENT:
            {
                // Fixed gradient: adjust owner cell based on gradient
                Scalar grad_n = bc->scalarValue();  // Normal gradient
                Scalar dist = face.d_Pf_mag();      // Cell to face distance

                // Face value with gradient: φ_face = φ_cell + grad_n * dist
                // Rearranging: φ_cell = φ_face - grad_n * dist
                // For now, apply gradient correction
                field[ownerCell] -= grad_n * dist;
                break;
            }

            case BCType::NO_SLIP:
            {
                // NO_SLIP is for velocity only, ignore for turbulence
                break;
            }

            default:
            {
                // Unknown BC type - use zero gradient
                break;
            }
        }
    }

    // Additional wall treatment for omega (viscous sublayer formula)
    // This overrides the BC values for cells very close to walls
    if (fieldName == "omega")
    {
        for (const auto& face : allFaces_)
        {
            if (!face.isBoundary()) continue;

            const BoundaryPatch* patch = faceToPatch.at(face.idx());

            // Only apply viscous sublayer formula at wall patches
            if (patch->type() == BoundaryConditionType::WALL)
            {
                size_t ownerCell = face.ownerCell();
                Scalar y = wallDistance_[ownerCell];

                // For viscous sublayer, use analytical formula
                // ω_wall = 6ν/(β₁y²)
                if (y < 1e-5)  // Very close to wall
                {
                    y = std::max(y, 1e-10);  // Prevent division by zero

                    field[ownerCell] = 
                        6.0 * nu_lam / (constants_.beta_1 * y * y);

                    // Limit maximum
                    field[ownerCell] = std::min(field[ownerCell], 1e7);
                }
            }
        }
    }

    // Wall treatment for k: k = 0 at wall face, but NOT at the owner cell!
    // The fixedValue BC from setup file already handles this correctly.
    // DO NOT override owner cell values here - that destroys the solution!
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
    
    Scalar y = wallDistance_[cellIdx];
    Scalar tau_wall = wallShearStress_[cellIdx];
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
        10.0 * constants_.beta_star * k_[cellIdx] * omega_[cellIdx];

    return std::min(P_k, limit);
}