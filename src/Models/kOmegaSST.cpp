/******************************************************************************
 * @file kOmegaSST.cpp
 * @brief Implementation of k-omega SST turbulence model
 *****************************************************************************/

#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

#include "kOmegaSST.h"
#include "Logger.h"
#include "Matrix.h"

// ************************* Constructor & Destructor *************************

kOmegaSST::kOmegaSST
(
    const Mesh& mesh,
    const BoundaryConditions& bc,
    const GradientScheme& gradientScheme,
    const ConvectionSchemes& kScheme,
    const ConvectionSchemes& omegaScheme
)
:
    mesh_(mesh),
    bcManager_(bc),
    gradientScheme_(gradientScheme),
    kConvectionScheme_(kScheme),
    omegaConvectionScheme_(omegaScheme),
    matrixConstruct_(std::make_unique<Matrix>(mesh_, bcManager_))
{
    if (debug_)
    {
        std::cout
            << "k-omega SST turbulence model initialized with "
            << mesh_.numCells() << " cells." << std::endl;
    }
}

kOmegaSST::~kOmegaSST() noexcept = default;


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

    yPlusLam();

    updateWallDistance();

    wallFunctionWeights();

    k_.setAll(initialK);
    omega_.setAll(initialOmega);

    boundTurbulenceFields();

    updateOmegaWallFunctionBoundaryValues();

    // Turbulent viscosity: nut = k / omega
    size_t numCells = mesh_.numCells();
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
    const TensorField& gradU
)
{
    boundTurbulenceFields();

    updateOmegaWallFunctionBoundaryValues();
    applyOmegaWallCellValues();

    const ScalarField strainRateField = strainRate(gradU);
    const VectorField gK = gradientK();
    const VectorField gO = gradientOmega();
    const ScalarField CDkOmega = crossDiffusion(gK, gO);

    const ScalarField f1 = F1(CDkOmega);
    const ScalarField f2 = F2();
    const ScalarField f3 = F3();
    const ScalarField f23 = F23(f2, f3);

    ProductionTerms P = productionTerms(f1, f23, strainRateField);
    overrideWallCellProduction(U, P.k);

    const ScalarField divUField = divU(flowRateFace);

    solveOmegaEquation(flowRateFace, f1, P.omega, CDkOmega, divUField, gO);
    boundOmega();

    solveKEquation(flowRateFace, f1, P.k, divUField, gK);
    boundK();
    boundOmega();

    updateTurbulentViscosity(f23, strainRateField);
    updateNutWallFace();
    updateWallShearStress(U);
    updateYPlus();

    logFieldDiagnostics();
}


// ****************************** Private Methods *****************************

void kOmegaSST::yPlusLam()
{
    Scalar yPlusLam = S(11.0);

    for (int i = 0; i < 10; ++i)
    {
        yPlusLam =
            std::log(std::max(coeffs_.E * yPlusLam, S(1.0)))
          / coeffs_.kappa;
    }

    yPlusLam_ = yPlusLam;
}

void kOmegaSST::updateWallDistance()
{
    if (debug_)
    {
        std::cout
            << "  Computing wall distance using meshWave method..."
            << std::endl;
    }

    // Initialize: large sentinel distance, zero wall point
    wallDistance_.setAll(S(1e10));
    nearestWallPoint_.setAll(Vector{});

    // Seed wall-adjacent cells with the perpendicular distance to each
    // wall face centroid
    for (const auto& face : mesh_.faces())
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch& patch = face.patch()->get();

        if (patch.type() != PatchType::WALL) continue;

        size_t cellIdx = face.ownerCell();
        Vector cellCenter = mesh_.cells()[cellIdx].centroid();
        Vector faceCenter = face.centroid();
        Vector normal = face.normal();

        Vector cellToFace = faceCenter - cellCenter;
        Scalar dist = std::abs(dot(cellToFace, normal));

        if (dist < wallDistance_[cellIdx])
        {
            wallDistance_[cellIdx] = dist;
            nearestWallPoint_[cellIdx] = faceCenter;
        }
    }

    // Iterative propagation: carry wall-point coordinates through
    // internal faces so each cell computes its own Euclidean distance
    const size_t maxIterations = 100;
    const Scalar tolerance = S(1e-12);

    for (size_t iter = 0; iter < maxIterations; ++iter)
    {
        Scalar maxChange = S(0.0);

        for (const auto& face : mesh_.faces())
        {
            if (face.isBoundary()) continue;

            size_t owner = face.ownerCell();
            size_t neighbor = face.neighborCell().value();

            Vector ownerCenter = mesh_.cells()[owner].centroid();
            Vector neighborCenter = mesh_.cells()[neighbor].centroid();

            // Try to improve owner using neighbor's nearest wall point
            {
                Vector candidatePoint = nearestWallPoint_[neighbor];
                Scalar newDist =
                    (ownerCenter - candidatePoint).magnitude();

                if (newDist < wallDistance_[owner])
                {
                    Scalar change = wallDistance_[owner] - newDist;
                    maxChange = std::max(maxChange, change);
                    wallDistance_[owner]     = newDist;
                    nearestWallPoint_[owner] = candidatePoint;
                }
            }

            // Try to improve neighbor using owner's nearest wall point
            {
                Vector candidatePoint = nearestWallPoint_[owner];
                Scalar newDist =
                    (neighborCenter - candidatePoint).magnitude();

                if (newDist < wallDistance_[neighbor])
                {
                    Scalar change = wallDistance_[neighbor] - newDist;
                    maxChange = std::max(maxChange, change);
                    wallDistance_[neighbor] = newDist;
                    nearestWallPoint_[neighbor] = candidatePoint;
                }
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

void kOmegaSST::wallFunctionWeights()
{
    const size_t numCells = mesh_.numCells();

    std::vector<size_t> wallFaceCountPerCell(numCells, 0);
    wallFaceWeight_.setAll(S(0.0));

    // Identify wall-function faces and accumulate wall area
    std::vector<Scalar> totalWallArea(numCells, S(0.0));
    wallFunctionFaceIndices_.clear();

    for (size_t faceIdx = 0; faceIdx < mesh_.numFaces(); ++faceIdx)
    {
        const auto& face = mesh_.faces()[faceIdx];
        if (!face.isBoundary()) continue;

        const BoundaryPatch& patch = face.patch()->get();
        if (patch.type() != PatchType::WALL) continue;

        const BoundaryData& bc =
            bcManager_.fieldBC(patch.patchName(), "omega");

        if (bc.type() != BCType::OMEGA_WALL_FUNCTION) continue;

        wallFunctionFaceIndices_.push_back(faceIdx);

        size_t cellIdx = face.ownerCell();
        totalWallArea[cellIdx] += face.projectedArea();
        ++wallFaceCountPerCell[cellIdx];
    }

    // Compute per-face weight = faceArea / totalWallArea
    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
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
        for (size_t faceIdx : mesh_.cells()[cellIdx].faceIndices())
        {
            const auto& face = mesh_.faces()[faceIdx];
            if (!face.isBoundary()) continue;

            const auto& patch = face.patch();
            if (patch.has_value()
                && patch->get().type() == PatchType::WALL)
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
    const FaceFluxField& flowRateFace,
    const ScalarField& f1,
    const ScalarField& productionOmega,
    const ScalarField& CDkOmega,
    const ScalarField& divUField,
    const VectorField& gradOmega
)
{
    size_t numCells = mesh_.numCells();

    ScalarField GammaOmega;
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar sigmaOmega =
            blend(f1[cellIdx], coeffs_.sigmaOmega1, coeffs_.sigmaOmega2);

        GammaOmega[cellIdx] = nu_ + sigmaOmega * nut_[cellIdx];
    }

    // Construct transport equation
    ScalarField omegaSource;

    TransportEquation equationOmega
    {
        .fieldName  = "omega",
        .phi        = omega_,
        .flowRate   = std::cref(flowRateFace),
        .convScheme = std::cref(omegaConvectionScheme_),
        .Gamma      = std::cref(GammaOmega),
        .GammaFace  = std::nullopt,
        .source     = omegaSource,
        .gradPhi    = gradOmega,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_->buildMatrix(equationOmega);

    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    // Add source terms
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar cellVolume = mesh_.cells()[cellIdx].volume();

        // Production term: Pk * gamma / nut
        vectorB(static_cast<Eigen::Index>(cellIdx)) +=
            productionOmega[cellIdx] * cellVolume;

        // Destruction term: -β·ω² (implicit: β·ω on diagonal)
        const Scalar beta =
            blend(f1[cellIdx], coeffs_.beta1, coeffs_.beta2);
        matrixA.coeffRef
        (
            static_cast<Eigen::Index>(cellIdx),
            static_cast<Eigen::Index>(cellIdx)
        ) += beta * omega_[cellIdx] * cellVolume;

        // Cross-diffusion: linearization of (1-F1)*CDkOmega
        Scalar CDkOmegaLineared =
            (S(1.0) - f1[cellIdx]) * CDkOmega[cellIdx]
          / (omega_[cellIdx] + vSmallValue);

        if (CDkOmegaLineared < S(0.0))
        {
            matrixA.coeffRef
            (
                static_cast<Eigen::Index>(cellIdx),
                static_cast<Eigen::Index>(cellIdx)
            ) += -CDkOmegaLineared * cellVolume;
        }
        else
        {
            vectorB(static_cast<Eigen::Index>(cellIdx)) +=
                CDkOmegaLineared * omega_[cellIdx] * cellVolume;
        }

        // -(2/3)*gamma*divU SuSp term (continuity correction)
        Scalar gamma = blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2);
        Scalar suspOmega =
            (S(2.0) / S(3.0)) * gamma * divUField[cellIdx];

        matrixA.coeffRef
        (
            static_cast<Eigen::Index>(cellIdx),
            static_cast<Eigen::Index>(cellIdx)
        ) += std::max(suspOmega, S(0.0)) * cellVolume;

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

    if (debug_)
    {
        Logger::residualRow
        (
            "omega",
            "BiCGSTAB",
            omegaSolver_.lastIterations(),
            omegaSolver_.lastResidual()
        );
    }
}

void kOmegaSST::updateOmegaWallFunctionBoundaryValues()
{
    omegaWallFunctionFaceValues_.setAll
    (
        std::numeric_limits<Scalar>::quiet_NaN()
    );

    constexpr Scalar minWallDistance = S(1e-20);

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        size_t cellIdx = face.ownerCell();

        const Scalar y =
            std::max
            (
                std::abs(dot(face.dPf(), face.normal())),
                minWallDistance
            );

        // Viscous sublayer
        Scalar omegaVis = S(6.0) * nu_ / (coeffs_.beta1 * y * y);

        // Log-region
        Scalar sqrtK =
            std::sqrt(std::max(k_[cellIdx], S(0.0)));

        Scalar omegaLog = sqrtK / (Cmu25_ * coeffs_.kappa * y + vSmallValue);

        Scalar yPlusFace = Cmu25_ * y * sqrtK / nu_;

        Scalar omegaWallValue = (yPlusFace < yPlusLam_) ? omegaVis : omegaLog;

        omegaWallFunctionFaceValues_[face.idx()] =
            std::max(omegaWallValue, omegaMin_);
    }
}

void kOmegaSST::applyOmegaWallCellValues()
{
    const size_t numCells = mesh_.numCells();

    std::vector<Scalar> omegaAccum(numCells, S(0.0));

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
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
    const FaceFluxField& flowRateFace,
    const ScalarField& f1,
    const ScalarField& productionK,
    const ScalarField& divUField,
    const VectorField& gradK
)
{
    size_t numCells = mesh_.numCells();

    // Construct transport matrix for k
    ScalarField GammaK;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar sigmaK = blend(f1[cellIdx], coeffs_.sigmaK1, coeffs_.sigmaK2);

        GammaK[cellIdx] = nu_ + sigmaK * nut_[cellIdx];
    }

    ScalarField kSource;

    TransportEquation equationK
    {
        .fieldName  = "k",
        .phi        = k_,
        .flowRate   = std::cref(flowRateFace),
        .convScheme = std::cref(kConvectionScheme_),
        .Gamma      = std::cref(GammaK),
        .GammaFace  = std::nullopt,
        .source     = kSource,
        .gradPhi    = gradK,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_->buildMatrix(equationK);

    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    // Add k source terms
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar cellVolume = mesh_.cells()[cellIdx].volume();

        // Apply Pk limiter
        Scalar Pk =
            std::min
            (
                productionK[cellIdx],
                coeffs_.c1 * coeffs_.betaStar * k_[cellIdx] * omega_[cellIdx]
            );

        vectorB(static_cast<Eigen::Index>(cellIdx)) += Pk * cellVolume;

        // Destruction term: -β*·kω
        const Scalar destruction = coeffs_.betaStar * omega_[cellIdx];
        matrixA.coeffRef
        (
            static_cast<Eigen::Index>(cellIdx),
            static_cast<Eigen::Index>(cellIdx)
        ) += destruction * cellVolume;

        // -(2/3)*divU SuSp term (continuity correction)
        Scalar suspK = (S(2.0) / S(3.0)) * divUField[cellIdx];
        matrixA.coeffRef
        (
            static_cast<Eigen::Index>(cellIdx),
            static_cast<Eigen::Index>(cellIdx)
        ) += std::max(suspK, S(0.0)) * cellVolume;

        vectorB(static_cast<Eigen::Index>(cellIdx)) +=
            std::max(-suspK, S(0.0)) * k_[cellIdx] * cellVolume;
    }

    // Apply implicit under-relaxation (Patankar's method)
    matrixConstruct_->relax(alphaK_, k_);

    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    kSolution(k_.data(), static_cast<Eigen::Index>(numCells));

    kSolver_.solveWithBiCGSTAB(kSolution, matrixA, vectorB);

    if (debug_)
    {
        Logger::residualRow
        (
            "k",
            "BiCGSTAB",
            kSolver_.lastIterations(),
            kSolver_.lastResidual()
        );
    }
}

void kOmegaSST::updateNutWallFace()
{
    nutWall_.setAll(S(0.0));

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        size_t cellIdx = face.ownerCell();

        // Wall-normal distance
        const Scalar y =
            std::max
            (
                std::abs(dot(face.dPf(), face.normal())),
                minWallDist_
            );

        Scalar sqrtK =
            std::sqrt(std::max(k_[cellIdx], S(0.0)));

        Scalar yPlusFace = Cmu25_ * sqrtK * y / nu_;

        if (yPlusFace > yPlusLam_)
        {
            // Log layer: nutw = nu * (yPlus*kappa/ln(E*yPlus) - 1)
            Scalar nutw =
                nu_
              * (
                    yPlusFace * coeffs_.kappa
                  / std::log(std::max(coeffs_.E * yPlusFace, S(1.0)))
                  - S(1.0)
                );

            nutWall_[face.idx()] = std::max(nutw, S(0.0));
        }
        // Viscous sublayer: nutw = 0 (already initialized)
    }
}

void kOmegaSST::updateTurbulentViscosity
(
    const ScalarField& f23,
    const ScalarField& strainRateField
)
{
    // SST turbulent viscosity:
    // nut = a1*k / max(a1*omega, b1*F23*sqrt(S2))
    size_t numCells = mesh_.numCells();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar strainMag = strainRateField[cellIdx];

        Scalar denominator =
            std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23[cellIdx] * strainMag
            );

        nut_[cellIdx] =
            coeffs_.a1 * k_[cellIdx] / (denominator + vSmallValue);

        nut_[cellIdx] = std::max(nut_[cellIdx], S(0.0));
    }
}

void kOmegaSST::updateWallShearStress
(
    const VectorField& U
)
{
    // Non-wall faces have zero shear stress
    wallShearStress_.setAll(S(0.0));

    constexpr Scalar maxKinematicShearStress = S(1000.0);

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        size_t cellIdx = face.ownerCell();

        Scalar y =
            std::max
            (
                std::abs(dot(face.dPf(), face.normal())),
                minWallDist_
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
        Scalar yPlusFace = Cmu25_ * sqrtK * y / nu_;

        Scalar tau;
        if (yPlusFace < yPlusLam_)
        {
            // Viscous sublayer: tau/rho = nu * U_tan / y
            tau = nu_ * U_tan_mag / (y + minWallDist_);
        }
        else
        {
            // Log layer: tau/rho = uTau^2 = Cmu^0.5 * k
            tau = Cmu25_ * Cmu25_ * k_[cellIdx];
        }

        wallShearStress_[face.idx()] =
            std::min(tau, maxKinematicShearStress);
    }
}

ScalarField kOmegaSST::F1(const ScalarField& CDkOmega) const
{
    const size_t numCells = mesh_.numCells();
    ScalarField f1;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar y = std::max(wallDistance_[cellIdx], vSmallValue);
        Scalar sqrtK = std::sqrt(k_[cellIdx]);

        // Read-only safe omega for division
        Scalar omegaSafe = std::max(omega_[cellIdx], omegaMin_);

        Scalar arg11 =
            sqrtK / (coeffs_.betaStar * omegaSafe * y);
        Scalar arg12 =
            S(500.0) * nu_ / (omegaSafe * y * y);

        // Cross-diffusion for blending function (clipped to positive)
        Scalar CDkw = std::max(CDkOmega[cellIdx], S(1e-10));

        Scalar arg13 =
            S(4.0) * coeffs_.sigmaOmega2 * k_[cellIdx]
          / (CDkw * y * y);

        Scalar arg1 =
            std::min
            (
                std::min(std::max(arg11, arg12), arg13),
                S(10.0)
            );

        Scalar arg1Sq = arg1 * arg1;
        f1[cellIdx] = std::tanh(arg1Sq * arg1Sq);
    }

    return f1;
}

ScalarField kOmegaSST::F2() const
{
    const size_t numCells = mesh_.numCells();
    ScalarField f2;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar y = std::max(wallDistance_[cellIdx], vSmallValue);
        Scalar sqrtK = std::sqrt(k_[cellIdx]);
        Scalar omegaSafe = std::max(omega_[cellIdx], omegaMin_);

        Scalar arg12 = S(500.0) * nu_ / (omegaSafe * y * y);

        Scalar arg2 =
            std::min
            (
                std::max
                (
                    S(2.0) * sqrtK
                  / (coeffs_.betaStar * omegaSafe * y),
                    arg12
                ),
                S(100.0)
            );

        f2[cellIdx] = std::tanh(arg2 * arg2);
    }

    return f2;
}

kOmegaSST::ProductionTerms kOmegaSST::productionTerms
(
    const ScalarField& f1,
    const ScalarField& f23,
    const ScalarField& strainRateField
) const
{
    size_t numCells = mesh_.numCells();

    ProductionTerms P;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar strainMag = strainRateField[cellIdx];
        Scalar strainSq = strainMag * strainMag;

        // SST Pk by nut limiter
        Scalar PkbyNutLimit =
            (coeffs_.c1 / coeffs_.a1) * coeffs_.betaStar * omega_[cellIdx]
          * std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23[cellIdx] * strainMag
            );

        Scalar PkByNut = std::min(strainSq, PkbyNutLimit);

        // k production = nut * S² (unlimited)
        P.k[cellIdx] = nut_[cellIdx] * strainSq;

        // omega production = gamma * GbyNu
        Scalar gamma =
            blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2);

        P.omega[cellIdx] = gamma * PkByNut;
    }

    return P;
}

void kOmegaSST::overrideWallCellProduction
(
    const VectorField& U,
    ScalarField& productionK
)
{
    const size_t numCells = mesh_.numCells();

    ScalarField wallProductionAccum;
    wallProductionAccum.setAll(S(0.0));

    std::vector<char> hasWallOverride(numCells, 0);

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
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
                minWallDist_
            );

        Scalar sqrtK =
            std::sqrt(std::max(k_[cellIdx], S(0.0)));

        Scalar yPlusFace = Cmu25_ * y * sqrtK / nu_;

        if (yPlusFace < yPlusLam_)
        {
            // Viscous sublayer contribute interior G so weights stay normalised
            wallProductionAccum[cellIdx] +=
                faceWeight * productionK[cellIdx];
            hasWallOverride[cellIdx] = 1;
            continue;
        }

        // G = sqr(uStar * magGradUw * y / uPlus) / (nu * kappa * yPlus)
        Scalar magGradUw = U[cellIdx].magnitude() / y;

        Scalar uStar = Cmu25_ * sqrtK;
        Scalar uPlus =
            std::log(std::max(coeffs_.E * yPlusFace, S(1.0)))
          / coeffs_.kappa;

        Scalar GwallArg =
            uStar * magGradUw * y / (uPlus + vSmallValue);

        Scalar Gwall =
            GwallArg * GwallArg
          / (nu_ * coeffs_.kappa * yPlusFace + vSmallValue);

        wallProductionAccum[cellIdx] += faceWeight * Gwall;
        hasWallOverride[cellIdx] = 1;
    }

    for (size_t i = 0; i < wallCellIndices_.size(); ++i)
    {
        size_t cellIdx = wallCellIndices_[i];
        if (!hasWallOverride[cellIdx]) continue;

        // Blend wall production with interior production using
        // wallCellFraction (wall area / total boundary area)
        Scalar f = wallCellFraction_[i];
        productionK[cellIdx] =
            std::lerp(productionK[cellIdx], wallProductionAccum[cellIdx], f);
    }
}

ScalarField kOmegaSST::strainRate
(
    const TensorField& gradU
) const
{
    size_t numCells = mesh_.numCells();
    ScalarField out;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // S = sqrt(2 * S_ij * S_ij) where S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        Scalar symmMagSq = gradU[cellIdx].symm().magnitudeSquared();
        out[cellIdx] = std::sqrt(S(2.0) * symmMagSq);
    }

    return out;
}

VectorField kOmegaSST::gradientK() const
{
    size_t numCells = mesh_.numCells();
    VectorField gK;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gK[cellIdx] = gradientScheme_.cellGradient("k", k_, cellIdx);
    }

    gradientScheme_.limitGradient("k", k_, gK);

    return gK;
}

VectorField kOmegaSST::gradientOmega() const
{
    size_t numCells = mesh_.numCells();
    VectorField gO;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gO[cellIdx] =
            gradientScheme_.cellGradient("omega", omega_, cellIdx);
    }

    gradientScheme_.limitGradient("omega", omega_, gO);

    return gO;
}

ScalarField kOmegaSST::crossDiffusion
(
    const VectorField& gK,
    const VectorField& gO
) const
{
    size_t numCells = mesh_.numCells();
    ScalarField CD;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        CD[cellIdx] =
            S(2.0) * coeffs_.sigmaOmega2
          * dot(gK[cellIdx], gO[cellIdx])
          / std::max(omega_[cellIdx], omegaMin_);
    }

    return CD;
}

ScalarField kOmegaSST::divU
(
    const FaceFluxField& flowRateFace
) const
{
    const size_t numCells = mesh_.numCells();
    ScalarField divUField;
    divUField.setAll(S(0.0));

    for (const auto& face : mesh_.faces())
    {
        Scalar flux = flowRateFace[face.idx()];
        size_t owner = face.ownerCell();

        if (face.isBoundary())
        {
            divUField[owner] += flux;
        }
        else
        {
            size_t neighbor = face.neighborCell().value();
            divUField[owner]    += flux;
            divUField[neighbor] -= flux;
        }
    }

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        divUField[cellIdx] /= mesh_.cells()[cellIdx].volume();
    }

    return divUField;
}

ScalarField kOmegaSST::F3() const
{
    const size_t numCells = mesh_.numCells();
    ScalarField f3;

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
        f3[cellIdx] = S(1.0) - std::tanh(arg3Sq * arg3Sq);
    }

    return f3;
}

ScalarField kOmegaSST::F23
(
    const ScalarField& f2,
    const ScalarField& f3
) const
{
    const size_t numCells = mesh_.numCells();
    ScalarField f23;

    if (useF3_)
    {
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            f23[cellIdx] = f2[cellIdx] * f3[cellIdx];
        }
    }
    else
    {
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            f23[cellIdx] = f2[cellIdx];
        }
    }

    return f23;
}

void kOmegaSST::updateYPlus()
{
    yPlus_.setAll(S(0.0));

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        size_t cellIdx = face.ownerCell();

        const Scalar y =
            std::max
            (
                std::abs(dot(face.dPf(), face.normal())),
                minWallDist_
            );

        Scalar sqrtK =
            std::sqrt(std::max(k_[cellIdx], S(0.0)));

        yPlus_[face.idx()] = Cmu25_ * sqrtK * y / nu_;
    }
}

void kOmegaSST::boundTurbulenceFields()
{
    const size_t numCells = mesh_.numCells();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        k_[cellIdx] = std::max(k_[cellIdx], kMin_);
        omega_[cellIdx] = std::max(omega_[cellIdx], omegaMin_);
    }
}

void kOmegaSST::boundOmega()
{
    const size_t numCells = mesh_.numCells();
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
    const size_t numCells = mesh_.numCells();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        k_[cellIdx] = std::max(k_[cellIdx], kMin_);
    }
}

void kOmegaSST::logFieldDiagnostics() const
{
    if (!debug_) return;

    size_t numCells = mesh_.numCells();
    if (numCells == 0) return;

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

    Logger::subsection("Turbulence field statistics");
    Logger::scalarStat("k",     kMin, kMax, kSum / n);
    Logger::scalarStat("omega", oMin, oMax, oSum / n);
    Logger::scalarStat("nut",   nMin, nMax, nSum / n);
}
