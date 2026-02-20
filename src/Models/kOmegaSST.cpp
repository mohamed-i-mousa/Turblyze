/******************************************************************************
 * @file kOmegaSST.cpp
 * @brief Implementation of k-omega SST turbulence model
 *****************************************************************************/

#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

#include "kOmegaSST.hpp"
#include "Matrix.hpp"
#include "ConvectionScheme.hpp"


// ************************* Constructor & Destructor *************************

kOmegaSST::kOmegaSST
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
      nut_("nut", cells.size(), 0.0),
      wallDistance_("wallDistance", cells.size(), 1.0),
      wallShearStress_("wallShearStress", cells.size(), 0.0),
      gradK_("gradK", cells.size(), Vector(0.0, 0.0, 0.0)),
      gradOmega_("gradOmega", cells.size(), Vector(0.0, 0.0, 0.0)),
      allFaces_(faces),
      allCells_(cells),
      bcManager_(bc),
      gradientScheme_(gradientScheme),
      kConvectionScheme_(kScheme),
      omegaConvectionScheme_(omegaScheme),
      F1_("F1", cells.size(), 0.0),
      F2_("F2", cells.size(), 0.0),
      productionK_("productionk", cells.size(), 0.0),
      productionOmega_("productionOmega", cells.size(), 0.0),
      crossDiffusion_("crossDiffusion", cells.size(), 0.0),
      strainRate_("strainRate", cells.size(), 0.0),
      omegaWallFunctionFaceValues_
      (
          "omegaWallFunctionFaceValues",
          faces.size(),
          std::numeric_limits<Scalar>::quiet_NaN()
      ),
      rho_(1.225),
      nu_(0.0),
      useF3_(false),
      kMin_(S(smallValue)),
      omegaMin_(S(smallValue)),
      alphaK_(0.5),
      alphaOmega_(0.5),
      kSolver_("k", S(1e-8), 1000),
      omegaSolver_("omega", S(1e-8), 1000)
{
    matrixConstruct_ =
        std::make_unique<Matrix>(allFaces_, allCells_, bcManager_);

    if (debug_)
    {
        std::cout
            << "k-omega SST turbulence model initialized with "
            << allCells_.size() << " cells." << std::endl;
    }
}

kOmegaSST::~kOmegaSST() = default;


// ****************************** Public Methods ******************************

void kOmegaSST::initialize
(
    Scalar nu,
    Scalar initialK,
    Scalar initialOmega,
    Scalar alphaK,
    Scalar alphaOmega
)
{
    nu_ = nu;
    alphaK_ = alphaK;
    alphaOmega_ = alphaOmega;

    if (debug_)
    {
        std::cout
            << std::endl
            << "=== Initializing k-omega SST Turbulence Model ===" 
            << std::endl;
    }

    calculateWallDistance();
    buildOmegaWallCellWeights();

    k_.setAll(initialK);
    omega_.setAll(initialOmega);

    boundTurbulenceFields();
    updateOmegaWallFunctionBoundaryValues();

    // Derive turbulent viscosity: nut = k / omega
    size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nut_[cellIdx] = k_[cellIdx] / (omega_[cellIdx] + vSmallValue);
        nut_[cellIdx] = std::max(nut_[cellIdx], S(0.0));
    }

    if (debug_)
    {
        std::cout
            << "Turbulence fields initialized successfully." << std::endl;
    }
}

void kOmegaSST::solve
(
    const VectorField& U,
    const FaceFluxField& flowRateFace,
    const std::vector<VectorField>& gradU
)
{
    if (omegaWallCellWeight_.size() != allCells_.size())
    {
        buildOmegaWallCellWeights();
    }

    boundTurbulenceFields();
    updateOmegaWallFunctionBoundaryValues();

    // computeGradK = true, computeGradOmega = true 
    calculateGradients(true, true);

    calculateBlendingFunctions();

    calculateStrainRate(gradU);

    calculateProductionTerms();

    solveOmegaEquation(flowRateFace);

    boundTurbulenceFields();

    updateOmegaWallFunctionBoundaryValues();

    // computeGradK = false, computeGradOmega = true 
    calculateGradients(false, true);

    calculateBlendingFunctions();

    calculateProductionTerms();

    overrideWallCellProduction(U);

    solveKEquation(flowRateFace);

    boundTurbulenceFields();

    calculateTurbulentViscosity();

    calculateWallShearStress(U);

    // Log turbulence field evolution diagnostics
    logFieldDiagnostics();
}

ScalarField kOmegaSST::getEffectiveViscosity() const
{
    size_t numCells = allCells_.size();
    ScalarField nuEff("nuEff", numCells);

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nuEff[cellIdx] = nu_ + nut_[cellIdx];
    }

    return nuEff;
}


// ****************************** Private Methods ******************************

void kOmegaSST::calculateWallDistance()
{
    if (debug_)
    {
        std::cout
            << "  Computing wall distance using meshWave method..."
            << std::endl;
    }

    // Initialize with large value
    wallDistance_.setAll(S(1e10));

    // Wall-adjacent cells
    for (const auto& face : allFaces_)
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch* patch = face.patch();

        if (patch->type() != BoundaryConditionType::WALL) continue;

        size_t cellIdx = face.ownerCell();
        Vector cellCenter = allCells_[cellIdx].centroid();
        Vector faceCenter = face.centroid();
        Vector normal = face.normal();

        Vector cellToFace = faceCenter - cellCenter;
        Scalar dist = std::abs(dot(cellToFace, normal));

        // Keep minimum distance if cell touches multiple walls
        wallDistance_[cellIdx] = std::min(wallDistance_[cellIdx], dist);
    }

    // Iterative propagation through internal faces
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

        if (maxChange < tolerance)
        {
            break;
        }

        if (iter == maxIterations - 1)
        {
            std::cout
                << "  WARNING: Reached max iterations without convergence"
                << std::endl;
        }
    }

    if (debug_)
    {
        std::cout
            << "  Wall distance calculation completed." << std::endl;
    }
}

void kOmegaSST::buildOmegaWallCellWeights()
{
    const size_t numCells = allCells_.size();

    omegaWallFaceCountPerCell_.assign(numCells, 0);
    omegaWallCellWeight_.assign(numCells, S(0.0));
    omegaMultiWallCellCount_ = 0;

    for (const auto& face : allFaces_)
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch* patch = face.patch();
        if (patch->type() != BoundaryConditionType::WALL)
        {
            continue;
        }

        const BoundaryData* bc =
            bcManager_.fieldBC(patch->patchName(), "omega");

        if (!bc || bc->type() != BCType::OMEGA_WALL_FUNCTION)
        {
            continue;
        }

        ++omegaWallFaceCountPerCell_[face.ownerCell()];
    }

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const size_t nFaces = omegaWallFaceCountPerCell_[cellIdx];
        if (nFaces == 0) continue;

        omegaWallCellWeight_[cellIdx] = S(1.0) / S(nFaces);

        if (nFaces > 1)
        {
            ++omegaMultiWallCellCount_;
        }
    }
}

void kOmegaSST::solveOmegaEquation
(
    const FaceFluxField& flowRateFace
)
{
    if (debug_)
    {
        std::cout
            << "  Solving omega transport equation (relaxation = "
            << alphaOmega_ << ")..." << std::endl;
    }

    size_t numCells = allCells_.size();

    ScalarField GammaOmega("GammaOmega", numCells);
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar sigmaOmega =
            blend(F1_[cellIdx], constants_.sigmaOmega1, constants_.sigmaOmega2);

        GammaOmega[cellIdx] = nu_ + sigmaOmega * nut_[cellIdx];
    }

    // Construct transport equation
    ScalarField omegaSource("omega_source", numCells, 0.0);

    TransportEquation equationOmega
    {
        omega_,                                    // phi
        "omega",                                   // fieldName
        std::cref(flowRateFace),                   // flowRate
        std::cref(omegaConvectionScheme_),         // convScheme
        std::cref(GammaOmega),                     // Gamma
        std::nullopt,                              // GammaFace
        omegaSource,                               // source
        gradOmega_,                                // gradPhi
        gradientScheme_,                           // gradScheme
        std::cref(omegaWallFunctionFaceValues_)    // boundaryFaceValues
    };

    matrixConstruct_->buildMatrix(equationOmega);

    auto& matrixA = matrixConstruct_->getMatrixA();
    auto& vectorB = matrixConstruct_->getVectorB();

    // Add source terms and modify diffusion
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar cellVolume = allCells_[cellIdx].volume();

        const Scalar beta = 
            blend(F1_[cellIdx], constants_.beta1, constants_.beta2);

        // Production: productionOmega_ already contains gamma * S^2
        // from calculateProductionTerms()
        Scalar POmega = productionOmega_[cellIdx];

        vectorB(cellIdx) += POmega * cellVolume;

        // Destruction term: -β·ω² (implicit: β·ω on diagonal)
        matrixA.coeffRef(cellIdx, cellIdx) +=
            beta * omega_[cellIdx] * cellVolume;

        // Cross-diffusion linearization -SuSp((F1-1)*CDkOmega/omega, omega)
        // mapped to this solver's A*phi=b convention:
        //   if susp > 0 : implicit sink      -> +susp*omega on diagonal
        //   if susp < 0 : explicit source    -> +(-susp*omega) on RHS
        Scalar dotKOmega = dot(gradK_[cellIdx], gradOmega_[cellIdx]);

        Scalar CDkOmega =
            std::max(2.0 * constants_.sigmaOmega2 * dotKOmega, S(1e-10));

        Scalar susp =
            (F1_[cellIdx] - S(1.0)) * CDkOmega
          / (omega_[cellIdx] + vSmallValue);

        const Scalar spCoeff = std::max(susp, S(0.0));

        const Scalar suCoeff = -std::min(susp, S(0.0));

        if (spCoeff > S(0.0))
        {
            matrixA.coeffRef(cellIdx, cellIdx) += spCoeff * cellVolume;
        }

        vectorB(cellIdx) +=
            suCoeff * omega_[cellIdx] * cellVolume;
    }

    // Apply implicit under-relaxation
    matrixConstruct_->relax(alphaOmega_, omega_);

    // Map omega field storage as Eigen vector (zero-copy)
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    omega_solution(omega_.data(), numCells);

    omegaSolver_.solveWithBiCGSTAB(omega_solution, matrixA, vectorB);

    // Clamp to minimum value
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        omega_[cellIdx] = std::max(omega_[cellIdx], omegaMin_);
    }

    // Wall function overrides near-wall cells with corner-weighted accumulation.
    std::vector<Scalar> omegaWallAccum(numCells, S(0.0));
    std::vector<char> hasOmegaWallValue(numCells, 0);

    for (const auto& face : allFaces_)
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch* patch = face.patch();
        if (patch->type() != BoundaryConditionType::WALL)
        {
            continue;
        }

        const BoundaryData* bc =
            bcManager_.fieldBC(patch->patchName(), "omega");

        if (!bc || bc->type() != BCType::OMEGA_WALL_FUNCTION)
        {
            continue;
        }

        const size_t cellIdx = face.ownerCell();
        const Scalar cellWeight = omegaWallCellWeight_[cellIdx];

        if (cellWeight <= S(0.0))
        {
            continue;
        }

        Scalar omegaWF = omegaWallFunctionFaceValues_[face.idx()];

        if (std::isfinite(omegaWF))
        {
            omegaWallAccum[cellIdx] += cellWeight * omegaWF;
            hasOmegaWallValue[cellIdx] = 1;
        }
        else if (!std::isnan(omegaWF))
        {
            std::cout
                << "  WARNING: Non-finite omega wall-function"
                << " value at face " << face.idx()
                << " (value = " << omegaWF << ")"
                << std::endl;
        }
    }

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        if (hasOmegaWallValue[cellIdx])
        {
            omega_[cellIdx] = std::max(omegaWallAccum[cellIdx], omegaMin_);
        }
    }

}

void kOmegaSST::updateOmegaWallFunctionBoundaryValues()
{
    if (debug_)
    {
        std::cout
            << "  Updating omega wall-function boundary values..."
            << std::endl;
    }

    omegaWallFunctionFaceValues_.setAll
    (
        std::numeric_limits<Scalar>::quiet_NaN()
    );

    for (const auto& face : allFaces_)
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch* patch = face.patch();

        if (patch->type() != BoundaryConditionType::WALL)
        {
            continue;
        }

        const BoundaryData* bc =
            bcManager_.fieldBC(patch->patchName(), "omega");

        if (!bc || bc->type() != BCType::OMEGA_WALL_FUNCTION)
        {
            continue;
        }

        size_t cellIdx = face.ownerCell();
        // Use wall-normal distance, not |dPf|, to avoid severe underprediction
        // on skewed faces where tangential offset dominates dPf magnitude.
        const Scalar y =
            std::max(std::abs(dot(face.d_Pf(), face.normal())), S(1e-20));

        // Viscous sublayer contribution
        Scalar omegaVis = 6.0 * nu_ / (constants_.beta1 * y * y);

        // Log-layer contribution
        Scalar Cmu25 = std::pow(constants_.Cmu, 0.25);
        Scalar omegaLog =
            std::sqrt(std::max(k_[cellIdx], S(0.0)))
          / (Cmu25 * constants_.kappa * y + vSmallValue);

        // Binomial2 blending
        Scalar omegaWallValue =
            std::sqrt(omegaVis * omegaVis + omegaLog * omegaLog);

        omegaWallFunctionFaceValues_[face.idx()] =
            std::max(omegaWallValue, omegaMin_);
    }
}

void kOmegaSST::solveKEquation
(
    const FaceFluxField& flowRateFace
)
{
    if (debug_)
    {
        std::cout
            << "  Solving k transport equation"
            << " (relaxation = "
            << alphaK_ << ")..." << std::endl;
    }

    size_t numCells = allCells_.size();

    // Construct transport matrix for k
    ScalarField GammaK("GammaK", numCells);

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar sigmaK =
            blend(F1_[cellIdx], constants_.sigmaK1, constants_.sigmaK2);

        GammaK[cellIdx] = nu_ + sigmaK * nut_[cellIdx];
    }

    ScalarField k_source("k_source", numCells, 0.0);

    TransportEquation equationK
    {
        k_,                            // phi
        "k",                           // fieldName
        std::cref(flowRateFace),       // flowRate
        std::cref(kConvectionScheme_), // convScheme
        std::cref(GammaK),            // Gamma
        std::nullopt,                  // GammaFace
        k_source,                      // source
        gradK_,                        // gradPhi
        gradientScheme_                // gradScheme
    };

    matrixConstruct_->buildMatrix(equationK);

    auto& A_matrix = matrixConstruct_->getMatrixA();
    auto& b_vector = matrixConstruct_->getVectorB();

    // Add k-specific source terms
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar cellVolume = allCells_[cellIdx].volume();

        // Production term from SST GbyNu limiter path
        b_vector(cellIdx) += productionK_[cellIdx] * cellVolume;

        // Destruction term: -β*·kω (kinematic formulation)
        const Scalar destruction = constants_.betaStar * omega_[cellIdx];
        A_matrix.coeffRef(cellIdx, cellIdx) += destruction * cellVolume;
    }

    // Apply implicit under-relaxation (Patankar's method)
    matrixConstruct_->relax(alphaK_, k_);

    // Map k field storage as Eigen vector (zero-copy)
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    k_solution(k_.data(), numCells);

    kSolver_.solveWithBiCGSTAB
    (
        k_solution,
        A_matrix,
        b_vector
    );

    // Clamp to minimum value
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        k_[cellIdx] = std::max(k_[cellIdx], kMin_);
    }

}

void kOmegaSST::calculateTurbulentViscosity()
{
    // SST turbulent viscosity:
    // nut = a1*k / max(a1*omega, b1*F23*sqrt(S2))
    if (debug_)
    {
        std::cout
            << "  Calculating turbulent viscosity..."
            << std::endl;
    }

    size_t numCells = allCells_.size();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar S = strainRate_[cellIdx];

        Scalar denominator =
            std::max
            (
                constants_.a1 * omega_[cellIdx],
                constants_.b1 * F23(cellIdx) * S
            );

        nut_[cellIdx] =
            constants_.a1 * k_[cellIdx] / (denominator + vSmallValue);

        nut_[cellIdx] = std::max(nut_[cellIdx], 0.0);
    }
}

void kOmegaSST::calculateWallShearStress
(
    const VectorField& U
)
{
    if (debug_)
    {
        std::cout
            << "  Calculating wall shear stress..." << std::endl;
    }

    // Non-wall cells have zero shear stress
    wallShearStress_.setAll(S(0.0));

    for (const auto& face : allFaces_)
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch* patch = face.patch();

        if (patch->type() != BoundaryConditionType::WALL)
        {
            continue;
        }

        size_t cellIdx = face.ownerCell();
        Scalar y = wallDistance_[cellIdx];

        // Project velocity onto wall plane (tangential component)
        Vector n = face.normal();
        Vector U_cell = U[cellIdx];
        Scalar U_normal = dot(U_cell, n);
        Vector U_tangential = U_cell - n * U_normal;
        Scalar U_tan_mag = U_tangential.magnitude();

        // Wall shear: tau = rho * nu * U_tangential / y
        wallShearStress_[cellIdx] =
            rho_ * nu_ * U_tan_mag / (y + vSmallValue);

        // Limit to reasonable values
        wallShearStress_[cellIdx] =
            std::min(wallShearStress_[cellIdx], S(1000.0));
    }
}

void kOmegaSST::calculateBlendingFunctions()
{
    size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar y = wallDistance_[cellIdx];
        Scalar sqrt_k = std::sqrt(k_[cellIdx]);

        // Avoid division by zero
        y = std::max(y, vSmallValue);
        omega_[cellIdx] = std::max(omega_[cellIdx], omegaMin_);

        // Arguments for blending functions
        Scalar arg11 = sqrt_k / (constants_.betaStar * omega_[cellIdx] * y);
        Scalar arg12 = 500.0 * nu_ / (omega_[cellIdx] * y * y);

        // Calculate CDkw for arg1_3 (cross-diffusion for blending function)
        Scalar dot_product = dot(gradK_[cellIdx], gradOmega_[cellIdx]);

        Scalar CDkw =
            std::max
            (
                2.0 * constants_.sigmaOmega2 / omega_[cellIdx] * dot_product,
                S(1e-10)
            );

        Scalar arg13 =
            4.0 * constants_.sigmaOmega2 * k_[cellIdx] / (CDkw * y * y);

        Scalar arg1 =
            std::min(std::min(std::max(arg11, arg12), arg13), 10.0);

        F1_[cellIdx] = std::tanh(std::pow(arg1, 4.0));

        Scalar arg2 = std::min
            (
                std::max
                (
                    2.0 * sqrt_k
                  / (constants_.betaStar * omega_[cellIdx] * y),
                    arg12
                ),
                100.0
            );

        F2_[cellIdx] = std::tanh(arg2 * arg2);
    }
}

void kOmegaSST::calculateProductionTerms()
{
    size_t numCells = allCells_.size();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar S = strainRate_[cellIdx];
        Scalar S2 = S * S;

        // SST GbyNu limiter
        Scalar GbyNuLimit =
            (constants_.c1 / constants_.a1)
          * constants_.betaStar * omega_[cellIdx]
          * std::max
            (
                constants_.a1 * omega_[cellIdx],
                constants_.b1 * F23(cellIdx) * S
            );

        Scalar GbyNu = std::min(S2, GbyNuLimit);

        // Production of k: P_k = nut * GbyNu
        productionK_[cellIdx] = nut_[cellIdx] * GbyNu;

        // Production of omega: P_omega = gamma * GbyNu
        Scalar gamma =
            blend(F1_[cellIdx], constants_.gamma1, constants_.gamma2);

        productionOmega_[cellIdx] = gamma * GbyNu;
    }
}

void kOmegaSST::overrideWallCellProduction
(
    const VectorField& U
)
{
    const Scalar Cmu25 = std::pow(constants_.Cmu, 0.25);
    const size_t numCells = allCells_.size();

    std::vector<Scalar> GwallAccum(numCells, S(0.0));
    std::vector<char> hasWallOverride(numCells, 0);

    for (const auto& face : allFaces_)
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch* patch = face.patch();
        if (patch->type() != BoundaryConditionType::WALL)
        {
            continue;
        }

        const BoundaryData* bc =
            bcManager_.fieldBC(patch->patchName(), "omega");

        if (!bc || bc->type() != BCType::OMEGA_WALL_FUNCTION)
        {
            continue;
        }

        size_t cellIdx = face.ownerCell();
        const Scalar cellWeight = omegaWallCellWeight_[cellIdx];

        if (cellWeight <= S(0.0))
        {
            continue;
        }

        // Wall-normal distance (same approach as omega wall function)
        const Scalar y =
            std::max
            (
                std::abs(dot(face.d_Pf(), face.normal())),
                S(1e-20)
            );

        // Tangential velocity at cell centre
        Vector n = face.normal();
        Vector Ucell = U[cellIdx];
        Vector Utan = Ucell - dot(Ucell, n) * n;
        Scalar magUtan = Utan.magnitude();

        // Wall velocity gradient
        Scalar magGradUw = magUtan / y;

        // Wall-function production of k (OpenFOAM v4 formula):
        // G = (nu + nut) * Cmu^0.25 * sqrt(k) * magGradUw
        //   / (kappa * y)
        Scalar Gwall =
            (nu_ + nut_[cellIdx])
          * Cmu25
          * std::sqrt(std::max(k_[cellIdx], S(0.0)))
          * magGradUw
          / (constants_.kappa * y);

        GwallAccum[cellIdx] += cellWeight * Gwall;
        hasWallOverride[cellIdx] = 1;
    }

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        if (!hasWallOverride[cellIdx]) continue;
        productionK_[cellIdx] = limitProduction(GwallAccum[cellIdx], cellIdx);
    }

}

void kOmegaSST::calculateStrainRate
(
    const std::vector<VectorField>& gradU
)
{
    size_t numCells = allCells_.size();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // S = sqrt(2 * S_ij * S_ij) where S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        Scalar S11 = gradU[0][cellIdx].x();
        Scalar S22 = gradU[1][cellIdx].y();
        Scalar S33 = gradU[2][cellIdx].z();
        Scalar S12 = 0.5 * (gradU[0][cellIdx].y() + gradU[1][cellIdx].x());
        Scalar S13 = 0.5 * (gradU[0][cellIdx].z() + gradU[2][cellIdx].x());
        Scalar S23 = 0.5 * (gradU[1][cellIdx].z() + gradU[2][cellIdx].y());

        strainRate_[cellIdx] = std::sqrt
            (
                2.0 * (S11*S11 + S22*S22 + S33*S33
              + 2.0*(S12*S12 + S13*S13 + S23*S23))
            );
    }
}

void kOmegaSST::calculateGradients
(
    bool computeGradK,
    bool computeGradOmega
)
{
    size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        if (computeGradK)
        {
            gradK_[cellIdx] = gradientScheme_.cellGradient("k",k_,cellIdx);
        }

        if (computeGradOmega)
        {
            gradOmega_[cellIdx] =
                gradientScheme_.cellGradient
                (
                    "omega",
                    omega_,
                    cellIdx,
                    &omegaWallFunctionFaceValues_
                );
        }
    }
}

void kOmegaSST::calculateCrossDiffusion()
{
    size_t numCells = allCells_.size();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar dot_product = dot(gradK_[cellIdx], gradOmega_[cellIdx]);

        crossDiffusion_[cellIdx] =
            2.0 * (1.0 - F1_[cellIdx]) * constants_.sigmaOmega2
          / (omega_[cellIdx] + vSmallValue) * dot_product;
    }
}

Scalar kOmegaSST::F3(size_t cellIdx) const
{
    Scalar y = std::max(wallDistance_[cellIdx], S(1e-12));
    Scalar omegaCell = std::max(omega_[cellIdx], omegaMin_);
    Scalar arg3 = 
        std::min(150.0 * nu_ / (omegaCell * y * y + vSmallValue), 10.0);
    return S(1.0) - std::tanh(std::pow(arg3, 4.0));
}

Scalar kOmegaSST::F23(size_t cellIdx) const
{
    if (!useF3_)
    {
        return F2_[cellIdx];
    }

    return F2_[cellIdx] * F3(cellIdx);
}

Scalar kOmegaSST::calculateYPlus(size_t cellIdx) const
{
    Scalar y = wallDistance_[cellIdx];
    Scalar tauWall = wallShearStress_[cellIdx];
    Scalar uTau = std::sqrt(tauWall / rho_);

    return uTau * y / nu_;
}

Scalar kOmegaSST::limitProduction(Scalar Pk, size_t cellIdx) const
{
    Scalar limit =
        10.0 * constants_.betaStar * k_[cellIdx] * omega_[cellIdx];

    return std::min(Pk, limit);
}

void kOmegaSST::boundTurbulenceFields()
{
    const size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        k_[cellIdx] = std::max(k_[cellIdx], kMin_);
        omega_[cellIdx] = std::max(omega_[cellIdx], omegaMin_);
    }
}

void kOmegaSST::logFieldDiagnostics() const
{
    if (!debug_) return;

    size_t numCells = allCells_.size();

    Scalar kMin = k_[0], kMax = k_[0], kSum = S(0.0);
    Scalar oMin = omega_[0], oMax = omega_[0], oSum = S(0.0);
    Scalar nMin = nut_[0], nMax = nut_[0], nSum = S(0.0);

    for (size_t i = 0; i < numCells; ++i)
    {
        kMin = std::min(kMin, k_[i]);
        kMax = std::max(kMax, k_[i]);
        kSum += k_[i];

        oMin = std::min(oMin, omega_[i]);
        oMax = std::max(oMax, omega_[i]);
        oSum += omega_[i];

        nMin = std::min(nMin, nut_[i]);
        nMax = std::max(nMax, nut_[i]);
        nSum += nut_[i];
    }

    Scalar n = S(numCells);

    std::cout
        << std::scientific
        << "  Turbulence: k     [min, max, mean] = ["
        << kMin << ", " << kMax << ", "
        << kSum / n << "]" << std::endl;

    std::cout
        << "  Turbulence: omega [min, max, mean] = ["
        << oMin << ", " << oMax << ", "
        << oSum / n << "]" << std::endl;

    std::cout
        << "  Turbulence: nut   [min, max, mean] = ["
        << nMin << ", " << nMax << ", "
        << nSum / n << "]" << std::endl;
}
