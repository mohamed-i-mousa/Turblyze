/******************************************************************************
 * @file kOmegaSST.cpp
 * @brief Implementation of k-omega SST turbulence model
 *****************************************************************************/

#include "kOmegaSST.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

#include <omp.h>

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
{}

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

    // Compute yPlusLam
    yPlusLam();

    // Compute wall distance field
    updateWallDistance();

    // Compute wall-function weights for boundary faces
    wallFunctionWeights();

    // Initialize turbulence fields with initial conditions
    k_.setAll(initialK);
    omega_.setAll(initialOmega);

    // initialize nut = k/omega
    initializeTurbulentViscosity();

    // Update y+ on wall faces
    updateYPlus();

    // Compute wall-function nut for first momentum solve
    updateNutWall();
}


void kOmegaSST::solve
(
    const VectorField& U,
    const FaceFluxField& flowRateFace,
    const TensorField& gradU
)
{
    // Update y+ on wall faces
    updateYPlus();

    // Compute goemtric quantites 
    const ScalarField strainRateField = strainRate(gradU);
    const ScalarField divUField = divU(flowRateFace);

    // Compute k Production
    kProduction(strainRateField);

    // Update wall-function boundary values for omega
    updateOmegaWallValues();

    // Pre-set wall-cell omega via area-weighted lerp
    applyOmegaWallCellValues();

    // Override k production at wall-adjacent cells
    overrideWallCellProduction(U, Pk_);

    // Compute gradients and cross-diffusion
    const VectorField gradK = gradientK();
    const VectorField gradOmega = gradientOmega();
    const ScalarField CDkOmega = crossDiffusion(gradK, gradOmega);

    // Compute blending functions
    const ScalarField f1 = F1(CDkOmega);
    const ScalarField f2 = F2();
    const ScalarField f3 = F3();
    const ScalarField f23 = F23(f2, f3);

    // Compute omega production
    omegaProduction(f1, strainRateField);

    // Apply SST production limiters
    limitProduction(Pk_, POmega_, f1, f23, strainRateField);

    // Solve omega transport equation
    solveOmegaEquation
    (
        flowRateFace,
        f1,
        POmega_,
        CDkOmega,
        divUField,
        gradOmega
    );
    boundOmega();

    // Solve k transport equation
    solveKEquation
    (
        flowRateFace,
        f1,
        Pk_,
        divUField,
        gradK
    );
    boundK();
    boundOmega();

    // Update turbulent viscosity with SST limiter
    updateTurbulentViscosity(f23, strainRateField);

    // Update wall-function nut on wall faces
    updateNutWall();

    // Update kinematic wall shear stress for diagnostics
    updateWallShearStress(U);

    // Log min/max/mean of k, omega, nut
    logFieldDiagnostics();
}

// ****************************** Private Methods *****************************

void kOmegaSST::yPlusLam()
{
    // Initial guess 
    Scalar yPlusLam = S(11.0);

    // 10 iterations to solve yPlusLam = E * yPlusLam / kappa * log(yPlusLam)
    for (int i = 0; i < 10; ++i)
    {
        yPlusLam =
            std::log(std::max(coeffs_.E * yPlusLam, S(1.0))) / coeffs_.kappa;
    }

    yPlusLam_ = yPlusLam;
}


void kOmegaSST::updateWallDistance()
{
    wallDistanceConverged_ = false;

    // Initialize: large sentinel distance, zero wall point
    wallDistance_.setAll(S(1e10));

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
            wallDistanceConverged_ = true;
            break;
        }
    }
}


void kOmegaSST::wallFunctionWeights()
{
    const size_t numCells = mesh_.numCells();

    std::vector<size_t> wallFaceCountPerCell(numCells, 0);
    std::vector<Scalar> totalWallArea(numCells, S(0.0));

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

        // cache the owner-cell wall-normal distance 
        y_[face.idx()] = std::max
        (
            std::abs(dot(face.dPf(), face.normal())),
            minWallDist_
        );
    }

    // Build unique wall cell indices
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

    #pragma omp parallel for schedule(static)
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
    const ScalarField& divU,
    const VectorField& gradOmega
)
{
    size_t numCells = mesh_.numCells();

    ScalarField GammaOmega;
    #pragma omp parallel for schedule(static)
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
    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar cellVolume = mesh_.cells()[cellIdx].volume();

        // Production term: Pk * gamma / nut
        vectorB(eIdx(cellIdx)) += productionOmega[cellIdx] * cellVolume;

        // Destruction term: -β·ω² (implicit: β·ω on diagonal)
        const Scalar beta = blend(f1[cellIdx], coeffs_.beta1, coeffs_.beta2);
        matrixA.coeffRef
        (
            eIdx(cellIdx),
            eIdx(cellIdx)
        ) += beta * omega_[cellIdx] * cellVolume;

        // Cross-diffusion: linearization of (1-F1)*CDkOmega
        Scalar CDkOmegaLineared =
            (S(1.0) - f1[cellIdx]) * CDkOmega[cellIdx]
          / (omega_[cellIdx] + vSmallValue);

        if (CDkOmegaLineared < S(0.0))
        {
            matrixA.coeffRef
            (
                eIdx(cellIdx),
                eIdx(cellIdx)
            ) += -CDkOmegaLineared * cellVolume;
        }
        else
        {
            vectorB(eIdx(cellIdx)) +=
                CDkOmegaLineared * omega_[cellIdx] * cellVolume;
        }

        // -(2/3)*gamma*divU SuSp term (continuity correction)
        Scalar gamma = blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2);
        Scalar suspOmega = (S(2.0) / S(3.0)) * gamma * divU[cellIdx];

        matrixA.coeffRef
        (
            eIdx(cellIdx),
            eIdx(cellIdx)
        ) += std::max(suspOmega, S(0.0)) * cellVolume;

        vectorB(eIdx(cellIdx)) +=
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
    omegaSolution(omega_.data(), eIdx(numCells));

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


void kOmegaSST::updateOmegaWallValues()
{
    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        size_t cellIdx = face.ownerCell();

        if (yPlus_[face.idx()] < yPlusLam_)
        {

            omegaWall_[face.idx()] =
                S(6.0) * nu_
              / (coeffs_.beta1 * y_[face.idx()] * y_[face.idx()]);
        }
        else
        {
            omegaWall_[face.idx()] = 
                std::sqrt(k_[cellIdx]) 
              / (Cmu25_ * coeffs_.kappa * y_[face.idx()]);
        }
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

        if (std::isfinite(omegaWall_[face.idx()]))
        {
            omegaAccum[cellIdx] += faceWeight * omegaWall_[face.idx()];
        }
    }

    for (size_t i = 0; i < wallCellIndices_.size(); ++i)
    {
        size_t cellIdx = wallCellIndices_[i];
        wallCellOmega_[i] = std::max(omegaAccum[cellIdx], omegaMin_);

        Scalar f = wallCellFraction_[i];
        omega_[cellIdx] = std::lerp(omega_[cellIdx], wallCellOmega_[i], f);
    }
}


void kOmegaSST::solveKEquation
(
    const FaceFluxField& flowRateFace,
    const ScalarField& f1,
    const ScalarField& productionK,
    const ScalarField& divU,
    const VectorField& gradK
)
{
    size_t numCells = mesh_.numCells();

    // Construct transport matrix for k
    ScalarField GammaK;

    #pragma omp parallel for schedule(static)
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
    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar cellVolume = mesh_.cells()[cellIdx].volume();

        vectorB(eIdx(cellIdx)) += productionK[cellIdx] * cellVolume;

        // Destruction term: -β*·kω
        const Scalar destruction = coeffs_.betaStar * omega_[cellIdx];
        matrixA.coeffRef(eIdx(cellIdx),eIdx(cellIdx)) +=
            destruction * cellVolume;

        // -(2/3)*divU SuSp term (continuity correction)
        Scalar suspK = (S(2.0) / S(3.0)) * divU[cellIdx];

        matrixA.coeffRef(eIdx(cellIdx),eIdx(cellIdx)) +=
            std::max(suspK, S(0.0)) * cellVolume;

        vectorB(eIdx(cellIdx)) +=
            std::max(-suspK, S(0.0)) * k_[cellIdx] * cellVolume;
    }

    // Apply implicit under-relaxation (Patankar's method)
    matrixConstruct_->relax(alphaK_, k_);

    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    kSolution(k_.data(), eIdx(numCells));

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


void kOmegaSST::updateNutWall()
{
    nutWall_.setAll(S(0.0));

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];

        if (yPlus_[face.idx()] > yPlusLam_)
        {
            // Log layer: nutw = nu * (yPlus*kappa/ln(E*yPlus) - 1)
            Scalar nutw =
                nu_
              * (
                    yPlus_[face.idx()] * coeffs_.kappa
                  / std::log(std::max(coeffs_.E * yPlus_[face.idx()], S(1.0)))
                  - S(1.0)
                );

            nutWall_[face.idx()] = std::max(nutw, S(0.0));
        }
        // Viscous sublayer: nutWall = 0 (already initialized)
    }
}


void kOmegaSST::initializeTurbulentViscosity()
{
    // Simple k-omega estimate: nut = k / omega
    // Used before strain rate & F23 are available for the SST limiter.
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nut_[cellIdx] =
            std::max(k_[cellIdx] / (omega_[cellIdx] + vSmallValue), S(0.0));
    }
}


void kOmegaSST::updateTurbulentViscosity
(
    const ScalarField& f23,
    const ScalarField& strainRate
)
{
    // SST turbulent viscosity:
    // nut = a1*k / max(a1*omega, b1*F23*sqrt(S2))
    size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nut_[cellIdx] =
            (coeffs_.a1 * k_[cellIdx])
          / std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23[cellIdx] * strainRate[cellIdx]
            );
    }
}


void kOmegaSST::updateWallShearStress
(
    const VectorField& U
)
{
    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        size_t cellIdx = face.ownerCell();

        // Project velocity onto wall plane (tangential component)
        Scalar normalVelocity = dot(U[cellIdx], face.normal());
        Vector tangentVelocity = U[cellIdx] - face.normal() * normalVelocity;
        Scalar tangentVelocityMag = tangentVelocity.magnitude();

        Scalar tau;
        if (yPlus_[face.idx()] < yPlusLam_)
        {
            // Viscous sublayer: tau/rho = nu * tangentVelocity / y
            tau = nu_ * tangentVelocityMag / y_[face.idx()];
        }
        else
        {
            // Log layer: tau/rho = uTau^2 = Cmu^0.5 * k
            tau = Cmu25_ * Cmu25_ * k_[cellIdx];
        }

        wallShearStress_[face.idx()] = std::min(tau, S(1000.0));
    }
}


ScalarField kOmegaSST::F1(const ScalarField& CDkOmega) const
{
    const size_t numCells = mesh_.numCells();
    ScalarField f1;

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar y = std::max(wallDistance_[cellIdx], vSmallValue);
        Scalar sqrtK = std::sqrt(k_[cellIdx]);

        // Cross-diffusion clipped to positive
        Scalar CDkw = std::max(CDkOmega[cellIdx], S(1e-10));

        Scalar arg1 =
            std::min
            (
                std::min
                (
                    std::max
                    (
                        sqrtK / (coeffs_.betaStar * omega_[cellIdx] * y),
                        S(500.0) * nu_ / (omega_[cellIdx] * y * y)
                    ),
                    S(4.0) * coeffs_.sigmaOmega2 * k_[cellIdx] / (CDkw * y * y)
                ),
                S(10.0)
            );

        // F1 = tanh(arg1^4)
        // Direct multiplication is faster that std::pow(arg1, 4)
        Scalar arg1Sq = arg1 * arg1;
        f1[cellIdx] = std::tanh(arg1Sq * arg1Sq);
    }

    return f1;
}


ScalarField kOmegaSST::F2() const
{
    const size_t numCells = mesh_.numCells();
    ScalarField f2;

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar y = std::max(wallDistance_[cellIdx], vSmallValue);
        Scalar sqrtK = std::sqrt(k_[cellIdx]);

        Scalar arg2 =
            std::min
            (
                std::max
                (
                    S(2.0) * sqrtK
                  / (coeffs_.betaStar * omega_[cellIdx] * y),
                    S(500.0) * nu_ / (omega_[cellIdx] * y * y)
                ),
                S(100.0)
            );

        // F2 = tanh(arg2^2)
        f2[cellIdx] = std::tanh(arg2 * arg2);
    }

    return f2;
}


ScalarField kOmegaSST::F3() const
{
    const size_t numCells = mesh_.numCells();
    ScalarField f3;

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar y = std::max(wallDistance_[cellIdx], S(1e-12));

        Scalar arg3 =
            std::min
            (
                S(150.0) * nu_ / (omega_[cellIdx] * y * y),
                S(10.0)
            );
         
        // F3 = 1 - tanh(arg3^4)    
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
        #pragma omp parallel for schedule(static)
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            f23[cellIdx] = f2[cellIdx] * f3[cellIdx];
        }
    }
    else
    {
        #pragma omp parallel for schedule(static)
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            f23[cellIdx] = f2[cellIdx];
        }
    }

    return f23;
}


void kOmegaSST::kProduction
(
    const ScalarField& strainRate
)
{
    size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar S2 = strainRate[cellIdx] * strainRate[cellIdx];

        // k production = nut * S² (unlimited)
        Pk_[cellIdx] = nut_[cellIdx] * S2;
    }
}


void kOmegaSST::omegaProduction
(
    const ScalarField& f1,
    const ScalarField& strainRate
)
{
    size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar S2 = strainRate[cellIdx] * strainRate[cellIdx];

        // omega production = gamma * GbyNut (unlimited)
        POmega_[cellIdx] = 
            blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2) * S2;
    }
}


void kOmegaSST::limitProduction
(
    ScalarField& productionK,
    ScalarField& productionOmega,
    const ScalarField& f1,
    const ScalarField& f23,
    const ScalarField& strainRate
)
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // Limit k production:
        Scalar kLimit =
            coeffs_.c1 * coeffs_.betaStar * k_[cellIdx] * omega_[cellIdx];
        productionK[cellIdx] = std::min(productionK[cellIdx], kLimit);

        // Limit omega production (Menter 2003 SST)
        Scalar omegaLimit =
            (coeffs_.c1 / coeffs_.a1) * coeffs_.betaStar * omega_[cellIdx]
          * blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2)
          * std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23[cellIdx] * strainRate[cellIdx]
            );
        productionOmega[cellIdx] =
            std::min(productionOmega[cellIdx], omegaLimit);
    }
}


void kOmegaSST::overrideWallCellProduction
(
    const VectorField& U,
    ScalarField& productionK
)
{
    const size_t numCells = mesh_.numCells();
    FaceData<Scalar> GWall;
    ScalarField wallProductionAccum;
    std::vector<char> hasWallOverride(numCells, 0);

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        size_t cellIdx = face.ownerCell();

        if (wallFaceWeight_[face.idx()] <= S(0.0))
        {
            continue;
        }

        if (yPlus_[face.idx()] < yPlusLam_)
        {
            // Viscous sublayer: contribute interior G
            GWall[face.idx()] = productionK[cellIdx];

            wallProductionAccum[cellIdx] +=
                wallFaceWeight_[face.idx()] * productionK[cellIdx];
            hasWallOverride[cellIdx] = 1;
        }
        else
        {
            // Log layer: G = sqr(uStar*magGradUw*y/uPlus) / (nu*kappa*yPlus)
            Scalar magGradUw = U[cellIdx].magnitude() / y_[face.idx()];
            Scalar uStar = Cmu25_ * std::sqrt(k_[cellIdx]);
            Scalar uPlus =
                std::log(std::max(coeffs_.E * yPlus_[face.idx()], S(1.0)))
              / coeffs_.kappa;
            Scalar uTau2 = uStar * magGradUw * y_[face.idx()] / (uPlus);

            GWall[face.idx()] =
                uTau2 * uTau2
              / (nu_ * coeffs_.kappa * yPlus_[face.idx()]);

            wallProductionAccum[cellIdx] +=
                wallFaceWeight_[face.idx()] * GWall[face.idx()];
            hasWallOverride[cellIdx] = 1;
        }
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
    ScalarField strainRateMag;

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // S = sqrt(2 * S_ij * S_ij) where S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        Scalar symmMagSq = gradU[cellIdx].symm().magnitudeSquared();
        strainRateMag[cellIdx] = std::sqrt(S(2.0) * symmMagSq);
    }

    return strainRateMag;
}

VectorField kOmegaSST::gradientK() const
{
    size_t numCells = mesh_.numCells();
    VectorField gradK;

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradK[cellIdx] = gradientScheme_.cellGradient("k", k_, cellIdx);
    }

    //gradientScheme_.limitGradient("k", k_, gradK);

    return gradK;
}


VectorField kOmegaSST::gradientOmega() const
{
    size_t numCells = mesh_.numCells();
    VectorField gradOmega;

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradOmega[cellIdx] =
            gradientScheme_.cellGradient("omega", omega_, cellIdx);
    }

    //gradientScheme_.limitGradient("omega", omega_, gradOmega);

    return gradOmega;
}


ScalarField kOmegaSST::crossDiffusion
(
    const VectorField& gradK,
    const VectorField& gradOmega
) const
{
    size_t numCells = mesh_.numCells();
    ScalarField CD;

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        CD[cellIdx] =
            S(2.0) * coeffs_.sigmaOmega2
          * dot(gradK[cellIdx], gradOmega[cellIdx])
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

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const auto& cell = mesh_.cells()[cellIdx];
        const auto& faceIndices = cell.faceIndices();
        const auto& faceSigns = cell.faceSigns();

        Scalar sum = S(0.0);
        for (size_t j = 0; j < faceIndices.size(); ++j)
        {
            sum += S(faceSigns[j]) * flowRateFace[faceIndices[j]];
        }

        divUField[cellIdx] = sum / cell.volume();
    }

    return divUField;
}


void kOmegaSST::updateYPlus()
{
    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        size_t cellIdx = face.ownerCell();

        yPlus_[face.idx()] =
            Cmu25_ * std::sqrt(k_[cellIdx]) * y_[face.idx()] / nu_;
    }
}


void kOmegaSST::boundOmega()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar omegaLowerBound =
            k_[cellIdx] / (maxViscosityRatio_ * nu_ + vSmallValue);

        omega_[cellIdx] =
            std::max(omega_[cellIdx], std::max(omegaMin_, omegaLowerBound));
    }
}


void kOmegaSST::boundK()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
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
