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

// ************************* Constructor & Destructor *************************

kOmegaSST::kOmegaSST
(
    std::span<const Face> faces,
    std::span<const Cell> cells,
    const BoundaryConditions& bc,
    const GradientScheme& gradientScheme,
    const ConvectionScheme& kScheme,
    const ConvectionScheme& omegaScheme
)
    : allFaces_(faces),
      allCells_(cells),
      bcManager_(bc),
      gradientScheme_(gradientScheme),
      kConvectionScheme_(kScheme),
      omegaConvectionScheme_(omegaScheme)
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

    calcYPlusLam();

    calculateWallDistance();

    buildWallFunctionWeights();

    k_.setAll(initialK);
    omega_.setAll(initialOmega);

    boundTurbulenceFields();

    updateOmegaWallFunctionBoundaryValues();

    // Turbulent viscosity: nut = k / omega
    size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nut_[cellIdx] = k_[cellIdx] / (omega_[cellIdx] + vSmallValue);
        nut_[cellIdx] = std::max(nut_[cellIdx], S(0.0));
    }

    // Compute wall-function nut for first momentum solve
    updateNutWallFace();

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
    std::span<const VectorField> gradU
)
{
    // 1. Bound k and omega
    boundTurbulenceFields();

    // 2. Wall treatment: set omega at wall faces and cells first
    updateOmegaWallFunctionBoundaryValues();
    applyOmegaWallCellValues();

    // 3. Strain rate from velocity gradients
    calculateStrainRate(gradU);

    // 4. Gradients, cross-diffusion, and blending functions
    calculateGradients(true, true);
    calculateCrossDiffusion();
    calculateBlendingFunctions();
    F3();
    F23();

    // 5. Production terms (use fresh F1/F2 and wall-updated omega)
    calculateProductionTerms();
    overrideWallCellProduction(U);

    // 6. Compute divU once for both equation solvers
    computeDivU(flowRateFace);

    // 7. Solve omega, then bound with viscosity-ratio limit
    solveOmegaEquation(flowRateFace);
    boundOmega();

    // 8. Solve k, then bound both k and omega
    solveKEquation(flowRateFace);
    boundK();
    boundOmega();

    // 9. Update turbulent viscosity and wall shear stress
    calculateTurbulentViscosity();
    updateNutWallFace();
    calculateWallShearStress(U);
    calculateYPlus();

    // Log turbulence field evolution diagnostics
    logFieldDiagnostics();
}

ScalarField kOmegaSST::effectiveViscosity() const
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

void kOmegaSST::calcYPlusLam()
{
    yPlusLam_ = S(11.0);

    for (int i = 0; i < 10; ++i)
    {
        yPlusLam_ =
            std::log(std::max(const_.E * yPlusLam_, S(1.0)))
          / const_.kappa;
    }
}

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

        if (patch->type() != PatchType::WALL) continue;

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
    const Scalar tolerance = S(1e-12);

    for (size_t iter = 0; iter < maxIterations; ++iter)
    {
        Scalar maxChange = S(0.0);

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

void kOmegaSST::buildWallFunctionWeights()
{
    const size_t numCells = allCells_.size();

    std::vector<size_t> wallFaceCountPerCell(numCells, 0);
    wallFaceWeight_.setAll(S(0.0));

    // Identify wall-function faces and accumulate wall area
    std::vector<Scalar> totalWallArea(numCells, S(0.0));
    wallFunctionFaceIndices_.clear();

    for (size_t faceIdx = 0; faceIdx < allFaces_.size(); ++faceIdx)
    {
        const auto& face = allFaces_[faceIdx];
        if (!face.isBoundary()) continue;

        const BoundaryPatch* patch = face.patch();
        if (patch->type() != PatchType::WALL) continue;

        const BoundaryData* bc =
            bcManager_.fieldBC(patch->patchName(), "omega");

        if (!bc || bc->type() != BCType::OMEGA_WALL_FUNCTION) continue;

        wallFunctionFaceIndices_.push_back(faceIdx);

        size_t cellIdx = face.ownerCell();
        totalWallArea[cellIdx] += face.projectedArea();
        ++wallFaceCountPerCell[cellIdx];
    }

    // Compute per-face weight = faceArea / totalWallArea
    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = allFaces_[faceIdx];
        size_t cellIdx = face.ownerCell();
        Scalar area = totalWallArea[cellIdx];

        if (area > S(0.0))
        {
            wallFaceWeight_[face.idx()] = face.projectedArea() / area;
        }
    }

    // Build unique wall cell indices
    wallCellIndices_.clear();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        if (wallFaceCountPerCell[cellIdx] > 0)
        {
            wallCellIndices_.push_back(cellIdx);
        }
    }

    // Compute wallCellFraction = wallFunctionArea / totalPolyWallArea
    std::vector<Scalar> totalPolyWallArea(numCells, S(0.0));

    for (size_t cellIdx : wallCellIndices_)
    {
        for (size_t faceIdx : allCells_[cellIdx].faceIndices())
        {
            const auto& face = allFaces_[faceIdx];
            if (!face.isBoundary()) continue;

            const BoundaryPatch* patch = face.patch();
            if (patch
                && patch->type() == PatchType::WALL)
            {
                totalPolyWallArea[cellIdx] += face.projectedArea();
            }
        }
    }

    wallCellFraction_.resize(wallCellIndices_.size());
    wallCellOmega_.assign(wallCellIndices_.size(), S(0.0));

    constexpr Scalar wallCellFractionTol = S(0.1);

    for (size_t i = 0; i < wallCellIndices_.size(); ++i)
    {
        size_t cellIdx = wallCellIndices_[i];

        Scalar rawFraction = S(1.0);
        if (totalPolyWallArea[cellIdx] > S(0.0))
        {
            rawFraction =
                totalWallArea[cellIdx] / totalPolyWallArea[cellIdx];
        }

        wallCellFraction_[i] =
            std::max
            (
                (rawFraction - wallCellFractionTol)
              / (S(1.0) - wallCellFractionTol),
                S(0.0)
            );
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
            blend(F1_[cellIdx], const_.sigmaOmega1, const_.sigmaOmega2);

        GammaOmega[cellIdx] = nu_ + sigmaOmega * nut_[cellIdx];
    }

    // Construct transport equation
    ScalarField omegaSource("omega_source", numCells, S(0.0));

    TransportEquation equationOmega
    {
        .fieldName  = "omega",
        .phi        = omega_,
        .flowRate   = std::cref(flowRateFace),
        .convScheme = std::cref(omegaConvectionScheme_),
        .Gamma      = std::cref(GammaOmega),
        .GammaFace  = std::nullopt,
        .source     = omegaSource,
        .gradPhi    = gradOmega_,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_->buildMatrix(equationOmega);

    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    // Add source terms
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar cellVolume = allCells_[cellIdx].volume();

        // Production term: Pk * gamma / nut
        Scalar POmega = productionOmega_[cellIdx];
        vectorB(static_cast<Eigen::Index>(cellIdx)) += POmega * cellVolume;

        // Destruction term: -β·ω² (implicit: β·ω on diagonal)
        const Scalar beta = blend(F1_[cellIdx], const_.beta1, const_.beta2);
        matrixA.coeffRef(static_cast<Eigen::Index>(cellIdx), static_cast<Eigen::Index>(cellIdx)) +=
            beta * omega_[cellIdx] * cellVolume;

        // Cross-diffusion: linearization of (1-F1)*CDkOmega
        Scalar CDkOmegaLineared =
            (S(1.0) - F1_[cellIdx]) * CDkOmega_[cellIdx]
          / (omega_[cellIdx] + vSmallValue);

        if (CDkOmegaLineared < S(0.0))
        {
            matrixA.coeffRef(static_cast<Eigen::Index>(cellIdx), static_cast<Eigen::Index>(cellIdx)) +=
                -CDkOmegaLineared * cellVolume;
        }
        else
        {
            vectorB(static_cast<Eigen::Index>(cellIdx)) +=
                CDkOmegaLineared * omega_[cellIdx] * cellVolume;
        }

        // -(2/3)*gamma*divU SuSp term (continuity correction)
        Scalar gamma = blend(F1_[cellIdx], const_.gamma1, const_.gamma2);
        Scalar suspOmega =
            (S(2.0) / S(3.0)) * gamma * divU_[cellIdx];

        matrixA.coeffRef(static_cast<Eigen::Index>(cellIdx), static_cast<Eigen::Index>(cellIdx)) +=
            std::max(suspOmega, S(0.0)) * cellVolume;

        vectorB(static_cast<Eigen::Index>(cellIdx)) +=
            std::max(-suspOmega, S(0.0)) * omega_[cellIdx] * cellVolume;
    }

    // Apply under-relaxation
    matrixConstruct_->relax(alphaOmega_, omega_);

    // Fix wall-cell rows to impose omega = omegaWall
    matrixConstruct_->setValues
    (
        wallCellIndices_,
        wallCellOmega_,
        wallCellFraction_
    );

    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    omegaSolution(omega_.data(), static_cast<Eigen::Index>(numCells));

    omegaSolver_.solveWithBiCGSTAB(omegaSolution, matrixA, vectorB);
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

    constexpr Scalar minWallDistance = S(1e-20);

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = allFaces_[faceIdx];
        size_t cellIdx = face.ownerCell();

        const Scalar y =
            std::max
            (
                std::abs(dot(face.dPf(), face.normal())),
                minWallDistance
            );

        // Viscous sublayer
        Scalar omegaVis = S(6.0) * nu_ / (const_.beta1 * y * y);

        // Log-region
        Scalar sqrtK =
            std::sqrt(std::max(k_[cellIdx], S(0.0)));

        Scalar omegaLog = sqrtK / (Cmu25_ * const_.kappa * y + vSmallValue);

        Scalar yPlus = Cmu25_ * y * sqrtK / nu_;

        Scalar omegaWallValue = (yPlus < yPlusLam_) ? omegaVis : omegaLog;

        omegaWallFunctionFaceValues_[face.idx()] =
            std::max(omegaWallValue, omegaMin_);
    }
}

void kOmegaSST::applyOmegaWallCellValues()
{
    const size_t numCells = allCells_.size();

    std::vector<Scalar> omegaAccum(numCells, S(0.0));

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = allFaces_[faceIdx];
        const size_t cellIdx = face.ownerCell();
        const Scalar faceWeight = wallFaceWeight_[face.idx()];

        if (faceWeight <= S(0.0)) continue;

        Scalar omegaWF = omegaWallFunctionFaceValues_[face.idx()];

        if (std::isfinite(omegaWF))
        {
            omegaAccum[cellIdx] += faceWeight * omegaWF;
        }
    }

    for (size_t i = 0; i < wallCellIndices_.size(); ++i)
    {
        size_t cellIdx = wallCellIndices_[i];
        wallCellOmega_[i] = std::max(omegaAccum[cellIdx], omegaMin_);

        Scalar f = wallCellFraction_[i];
        omega_[cellIdx] =
            std::lerp(omega_[cellIdx], wallCellOmega_[i], f);
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
        Scalar sigmaK = blend(F1_[cellIdx], const_.sigmaK1, const_.sigmaK2);

        GammaK[cellIdx] = nu_ + sigmaK * nut_[cellIdx];
    }

    ScalarField k_source("k_source", numCells, S(0.0));

    TransportEquation equationK
    {
        .fieldName  = "k",
        .phi        = k_,
        .flowRate   = std::cref(flowRateFace),
        .convScheme = std::cref(kConvectionScheme_),
        .Gamma      = std::cref(GammaK),
        .GammaFace  = std::nullopt,
        .source     = k_source,
        .gradPhi    = gradK_,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_->buildMatrix(equationK);

    auto& A_matrix = matrixConstruct_->matrixA();
    auto& b_vector = matrixConstruct_->vectorB();

    // Add k source terms
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar cellVolume = allCells_[cellIdx].volume();

        // Apply Pk limiter
        Scalar Pk =
            std::min
            (
                productionK_[cellIdx],
                const_.c1 * const_.betaStar * k_[cellIdx] * omega_[cellIdx]
            );

        b_vector(static_cast<Eigen::Index>(cellIdx)) += Pk * cellVolume;

        // Destruction term: -β*·kω
        const Scalar destruction = const_.betaStar * omega_[cellIdx];
        A_matrix.coeffRef(static_cast<Eigen::Index>(cellIdx), static_cast<Eigen::Index>(cellIdx)) += destruction * cellVolume;

        // -(2/3)*divU SuSp term (continuity correction)
        Scalar suspK = (S(2.0) / S(3.0)) * divU_[cellIdx];
        A_matrix.coeffRef(static_cast<Eigen::Index>(cellIdx), static_cast<Eigen::Index>(cellIdx)) +=
            std::max(suspK, S(0.0)) * cellVolume;

        b_vector(static_cast<Eigen::Index>(cellIdx)) +=
            std::max(-suspK, S(0.0)) * k_[cellIdx] * cellVolume;
    }

    // Apply implicit under-relaxation (Patankar's method)
    matrixConstruct_->relax(alphaK_, k_);

    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    k_solution(k_.data(), static_cast<Eigen::Index>(numCells));

    kSolver_.solveWithBiCGSTAB
    (
        k_solution,
        A_matrix,
        b_vector
    );
}

void kOmegaSST::updateNutWallFace()
{
    nutWall_.setAll(S(0.0));

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = allFaces_[faceIdx];
        size_t cellIdx = face.ownerCell();

        // Wall-normal distance
        const Scalar y =
            std::max
            (
                std::abs(dot(face.dPf(), face.normal())),
                vSmallValue
            );

        Scalar sqrtK =
            std::sqrt(std::max(k_[cellIdx], S(0.0)));

        Scalar yPlus = Cmu25_ * sqrtK * y / nu_;

        if (yPlus > yPlusLam_)
        {
            // Log layer: nutw = nu * (yPlus*kappa/ln(E*yPlus) - 1)
            Scalar nutw =
                nu_
              * (
                    yPlus * const_.kappa
                  / std::log(std::max(const_.E * yPlus, S(1.0)))
                  - S(1.0)
                );

            nutWall_[face.idx()] = std::max(nutw, S(0.0));
        }
        // Viscous sublayer: nutw = 0 (already initialized)
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
        Scalar strainMag = strainRate_[cellIdx];

        Scalar denominator =
            std::max
            (
                const_.a1 * omega_[cellIdx],
                const_.b1 * F23_[cellIdx] * strainMag
            );

        nut_[cellIdx] =
            const_.a1 * k_[cellIdx] / (denominator + vSmallValue);

        nut_[cellIdx] = std::max(nut_[cellIdx], S(0.0));
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

    // Non-wall faces have zero shear stress
    wallShearStress_.setAll(S(0.0));

    constexpr Scalar minWallDistance = S(1e-20);
    constexpr Scalar maxWallShearStress = S(1000.0);

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = allFaces_[faceIdx];
        size_t cellIdx = face.ownerCell();

        Scalar y =
            std::max
            (
                std::abs(dot(face.dPf(), face.normal())),
                minWallDistance
            );

        // Project velocity onto wall plane (tangential component)
        Vector n = face.normal();
        Vector U_cell = U[cellIdx];
        Scalar U_normal = dot(U_cell, n);
        Vector U_tangential = U_cell - n * U_normal;
        Scalar U_tan_mag = U_tangential.magnitude();

        // Compute y+ to determine layer treatment
        Scalar sqrtK =
            std::sqrt(std::max(k_[cellIdx], S(0.0)));
        Scalar yPlus = Cmu25_ * sqrtK * y / nu_;

        Scalar tau;
        if (yPlus < yPlusLam_)
        {
            // Viscous sublayer: tau = rho * nu * U_tan / y
            tau = rho_ * nu_ * U_tan_mag / (y + vSmallValue);
        }
        else
        {
            // Log layer: tau = rho * uTau^2 = rho * Cmu^0.5 * k
            tau = rho_ * Cmu25_ * Cmu25_ * k_[cellIdx];
        }

        wallShearStress_[face.idx()] =
            std::min(tau, maxWallShearStress);
    }
}

void kOmegaSST::calculateBlendingFunctions()
{
    size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar y = wallDistance_[cellIdx];
        Scalar sqrtK = std::sqrt(k_[cellIdx]);

        // Avoid division by zero
        y = std::max(y, vSmallValue);
        omega_[cellIdx] = std::max(omega_[cellIdx], omegaMin_);

        // Arguments for blending functions
        Scalar arg11 =
            sqrtK / (const_.betaStar * omega_[cellIdx] * y);
        Scalar arg12 =
            S(500.0) * nu_ / (omega_[cellIdx] * y * y);

        // Cross-diffusion for blending function (clipped to positive)
        Scalar CDkw = std::max(CDkOmega_[cellIdx], S(1e-10));

        Scalar arg13 =
            S(4.0) * const_.sigmaOmega2 * k_[cellIdx]
          / (CDkw * y * y);

        Scalar arg1 =
            std::min
            (
                std::min(std::max(arg11, arg12), arg13),
                S(10.0)
            );

        Scalar arg1Sq = arg1 * arg1;
        F1_[cellIdx] = std::tanh(arg1Sq * arg1Sq);

        Scalar arg2 = std::min
            (
                std::max
                (
                    S(2.0) * sqrtK
                  / (const_.betaStar * omega_[cellIdx] * y),
                    arg12
                ),
                S(100.0)
            );

        F2_[cellIdx] = std::tanh(arg2 * arg2);
    }
}

void kOmegaSST::calculateProductionTerms()
{
    size_t numCells = allCells_.size();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar strainMag = strainRate_[cellIdx];
        Scalar strainSq = strainMag * strainMag;

        // SST Pk by nut limiter
        Scalar PkbyNutLimit =
            (const_.c1 / const_.a1) * const_.betaStar * omega_[cellIdx]
          * std::max
            (
                const_.a1 * omega_[cellIdx],
                const_.b1 * F23_[cellIdx] * strainMag
            );

        Scalar PkByNut = std::min(strainSq, PkbyNutLimit);

        // k production = nut * S² (unlimited)
        productionK_[cellIdx] = nut_[cellIdx] * strainSq;

        // omega production = gamma * GbyNu
        Scalar gamma =
            blend(F1_[cellIdx], const_.gamma1, const_.gamma2);

        productionOmega_[cellIdx] = gamma * PkByNut;
    }
}

void kOmegaSST::overrideWallCellProduction
(
    const VectorField& U
)
{
    const size_t numCells = allCells_.size();

    std::vector<Scalar> GwallAccum(numCells, S(0.0));
    std::vector<char> hasWallOverride(numCells, 0);

    constexpr Scalar minWallDistance = S(1e-20);

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = allFaces_[faceIdx];
        size_t cellIdx = face.ownerCell();
        const Scalar faceWeight = wallFaceWeight_[face.idx()];

        if (faceWeight <= S(0.0))
        {
            continue;
        }

        const Scalar y =
            std::max
            (
                std::abs(dot(face.dPf(), face.normal())),
                minWallDistance
            );

        Scalar sqrtK =
            std::sqrt(std::max(k_[cellIdx], S(0.0)));

        Scalar yPlus = Cmu25_ * y * sqrtK / nu_;

        // Viscous sublayer: keep interior G unchanged
        if (yPlus < yPlusLam_)
        {
            continue;
        }

        // G = sqr(uStar * magGradUw * y / uPlus) / (nu * kappa * yPlus)
        Scalar magGradUw = U[cellIdx].magnitude() / y;

        Scalar uStar = Cmu25_ * sqrtK;
        Scalar uPlus =
            std::log(std::max(const_.E * yPlus, S(1.0)))
          / const_.kappa;

        Scalar GwallArg =
            uStar * magGradUw * y / (uPlus + vSmallValue);

        Scalar Gwall =
            GwallArg * GwallArg
          / (nu_ * const_.kappa * yPlus + vSmallValue);

        GwallAccum[cellIdx] += faceWeight * Gwall;
        hasWallOverride[cellIdx] = 1;
    }

    for (size_t i = 0; i < wallCellIndices_.size(); ++i)
    {
        size_t cellIdx = wallCellIndices_[i];
        if (!hasWallOverride[cellIdx]) continue;

        // Blend wall production with interior production using
        // wallCellFraction (wall area / total boundary area)
        Scalar f = wallCellFraction_[i];
        productionK_[cellIdx] =
            std::lerp(productionK_[cellIdx], GwallAccum[cellIdx], f);
    }
}

void kOmegaSST::calculateStrainRate
(
    std::span<const VectorField> gradU
)
{
    size_t numCells = allCells_.size();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // S = sqrt(2 * S_ij * S_ij) where S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        Scalar S11 = gradU[0][cellIdx].x();
        Scalar S22 = gradU[1][cellIdx].y();
        Scalar S33 = gradU[2][cellIdx].z();
        Scalar S12 = S(0.5) * (gradU[0][cellIdx].y() + gradU[1][cellIdx].x());
        Scalar S13 = S(0.5) * (gradU[0][cellIdx].z() + gradU[2][cellIdx].x());
        Scalar S23 = S(0.5) * (gradU[1][cellIdx].z() + gradU[2][cellIdx].y());

        strainRate_[cellIdx] = std::sqrt
            (
                S(2.0) * (S11*S11 + S22*S22 + S33*S33
              + S(2.0) * (S12*S12 + S13*S13 + S23*S23))
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
                    cellIdx
                );
        }
    }

    // Apply Barth-Jespersen cell-based gradient limiting
    if (computeGradK)
    {
        gradientScheme_.limitGradient("k", k_, gradK_);
    }
    if (computeGradOmega)
    {
        gradientScheme_.limitGradient("omega", omega_, gradOmega_);
    }
}

void kOmegaSST::calculateCrossDiffusion()
{
    size_t numCells = allCells_.size();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        CDkOmega_[cellIdx] =
            S(2.0) * const_.sigmaOmega2
          * dot(gradK_[cellIdx], gradOmega_[cellIdx])
          / (omega_[cellIdx] + vSmallValue);
    }
}

void kOmegaSST::computeDivU
(
    const FaceFluxField& flowRateFace
)
{
    const size_t numCells = allCells_.size();
    divU_.setAll(S(0.0));

    for (const auto& face : allFaces_)
    {
        Scalar flux = flowRateFace[face.idx()];
        size_t owner = face.ownerCell();

        if (face.isBoundary())
        {
            divU_[owner] += flux;
        }
        else
        {
            size_t neighbor = face.neighborCell().value();
            divU_[owner]    += flux;
            divU_[neighbor] -= flux;
        }
    }

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        divU_[cellIdx] /= allCells_[cellIdx].volume();
    }
}


void kOmegaSST::F3()
{
    const size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar y = std::max(wallDistance_[cellIdx], S(1e-12));
        Scalar omegaCell = std::max(omega_[cellIdx], omegaMin_);
        Scalar arg3 =
            std::min
            (
                S(150.0) * nu_ / (omegaCell * y * y + vSmallValue),
                S(10.0)
            );
        Scalar arg3Sq = arg3 * arg3;
        F3_[cellIdx] = S(1.0) - std::tanh(arg3Sq * arg3Sq);
    }
}

void kOmegaSST::F23()
{
    const size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        if (!useF3_)
        {
            F23_[cellIdx] = F2_[cellIdx];
        }
        else
        {
            F23_[cellIdx] = F2_[cellIdx] * F3_[cellIdx];
        }
    }
}

void kOmegaSST::calculateYPlus()
{
    yPlus_.setAll(S(0.0));

    constexpr Scalar minWallDistance = S(1e-20);

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = allFaces_[faceIdx];
        size_t cellIdx = face.ownerCell();

        const Scalar y =
            std::max
            (
                std::abs(dot(face.dPf(), face.normal())),
                minWallDistance
            );

        Scalar sqrtK =
            std::sqrt(std::max(k_[cellIdx], S(0.0)));

        yPlus_[face.idx()] = Cmu25_ * sqrtK * y / nu_;
    }
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

void kOmegaSST::boundOmega()
{
    const size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar omegaLowerBound =
            k_[cellIdx] / (nutMaxCoeff_ * nu_ + vSmallValue);

        omega_[cellIdx] =
            std::max(omega_[cellIdx], std::max(omegaMin_, omegaLowerBound));
    }
}

void kOmegaSST::boundK()
{
    const size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        k_[cellIdx] = std::max(k_[cellIdx], kMin_);
    }
}

void kOmegaSST::logFieldDiagnostics() const
{
    if (!debug_) return;

    size_t numCells = allCells_.size();

    Scalar kMin = k_[0];
    Scalar kMax = k_[0];
    Scalar kSum = S(0.0);

    Scalar oMin = omega_[0];
    Scalar oMax = omega_[0];
    Scalar oSum = S(0.0);
    
    Scalar nMin = nut_[0];
    Scalar nMax = nut_[0];
    Scalar nSum = S(0.0);

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
        << kSum / n << ']' << std::endl;

    std::cout
        << "  Turbulence: omega [min, max, mean] = ["
        << oMin << ", " << oMax << ", "
        << oSum / n << ']' << std::endl;

    std::cout
        << "  Turbulence: nut   [min, max, mean] = ["
        << nMin << ", " << nMax << ", "
        << nSum / n << ']' << std::endl;
}
