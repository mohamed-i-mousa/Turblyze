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

// ************************* Special Member Functions *************************

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

    hasWallOverride_.assign(mesh_.numCells(), 0);
    omegaAccum_.assign(mesh_.numCells(), S(0.0));

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
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz,
    const FaceFluxField& flowRateFace,
    const TensorField& gradU
)
{
    if (!kSolver_)
    {
        FatalError
        (
            "kOmegaSST::solve: k linear solver not set; "
            "call setKSolver() (or SIMPLE::setTurbulenceSolvers()) "
            "before solve()."
        );
    }
    if (!omegaSolver_)
    {
        FatalError
        (
            "kOmegaSST::solve: omega linear solver not set; "
            "call setOmegaSolver() (or SIMPLE::setTurbulenceSolvers()) "
            "before solve()."
        );
    }

    // Update y+ on wall faces
    updateYPlus();

    // Compute geometric quantities
    strainRate(gradU);
    divU(flowRateFace);

    // Compute k Production
    kProduction();

    // Update wall-function boundary values for omega
    updateOmegaWallValues();

    // Pre-set wall-cell omega via area-weighted lerp
    applyOmegaWallCellValues();

    // Override k production at wall-adjacent cells
    overrideWallCellProduction(Ux, Uy, Uz);

    // Compute gradients and cross-diffusion
    gradientScheme_.fieldGradient(Field::k, k_, gradK_);
    gradientScheme_.fieldGradient(Field::omega, omega_, gradOmega_);

    crossDiffusion();

    // Compute blending functions
    F1();
    F2();
    F3();
    F23();

    // Compute omega production
    omegaProduction();

    // Apply SST production limiters
    limitProduction();

    // Solve omega transport equation
    solveOmegaEquation(flowRateFace);
    boundOmega();

    // Solve k transport equation
    solveKEquation(flowRateFace);
    boundK();
    boundOmega();

    // Update turbulent viscosity with SST limiter
    updateTurbulentViscosity();

    // Update wall-function nut on wall faces
    updateNutWall();

    // Update kinematic wall shear stress for diagnostics
    updateWallShearStress(Ux, Uy, Uz);

    // Log min/max/mean of k, omega, nut
    logFieldDiagnostics();
}


Scalar kOmegaSST::inletK
(
    const Vector& velocity,
    Scalar turbulenceIntensity
) noexcept
{
    const Scalar uPrime = turbulenceIntensity * magnitude(velocity);
    return std::max(S(1.5) * uPrime * uPrime, smallValue);
}


Scalar kOmegaSST::inletOmega
(
    Scalar k,
    Scalar hydraulicDiameter
) noexcept
{
    const Scalar lengthScale = S(0.07) * hydraulicDiameter;

    const Scalar omegaValue =
        std::sqrt(std::max(k, S(0.0)))
      / (std::pow(coeffs_.betaStar, S(0.25)) * lengthScale);

    return std::max(omegaValue, smallValue);
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
    wallDistance_.setAll(S(1e10));
    nearestWallPoint_.setAll(Vector{S(1e15), S(1e15), S(1e15)});

    // Seed wall-adjacent cells with the perpendicular distance to each
    // wall face centroid
    for (const auto& face : mesh_.faces())
    {
        if (!face.isBoundary()) continue;

        const BoundaryPatch& patch = face.patch()->get();

        if (patch.type() != PatchType::wall) continue;

        const size_t cellIdx = face.ownerCell();
        const Vector cellCenter = mesh_.cells()[cellIdx].centroid();
        const Vector faceCenter = face.centroid();
        const Vector normal = face.normal();
        const Vector cellToFace = faceCenter - cellCenter;
        const Scalar dist = std::abs(dot(cellToFace, normal));

        if (dist < wallDistance_[cellIdx])
        {
            wallDistance_[cellIdx] = dist;
            nearestWallPoint_[cellIdx] = faceCenter;
        }
    }

    // Iterative propagation: carry wall-point coordinates through
    // internal faces so each cell computes its own Euclidean distance
    const size_t maxIterations = 100;
    const Scalar tolerance = smallValue;

    for (size_t iter = 0; iter < maxIterations; ++iter)
    {
        Scalar maxChange = S(0.0);

        for (const auto& face : mesh_.faces())
        {
            if (face.isBoundary()) continue;

            const size_t owner = face.ownerCell();
            const size_t neighbor = face.neighborCell().value();
            const Vector ownerCenter = mesh_.cells()[owner].centroid();
            const Vector neighborCenter = mesh_.cells()[neighbor].centroid();

            // Try to improve owner using neighbor's nearest wall point
            {
                const Vector candidatePoint = nearestWallPoint_[neighbor];
                const Scalar newDist =
                    magnitude(ownerCenter - candidatePoint);

                if (newDist < wallDistance_[owner])
                {
                    const Scalar change = wallDistance_[owner] - newDist;
                    maxChange = std::max(maxChange, change);
                    wallDistance_[owner]     = newDist;
                    nearestWallPoint_[owner] = candidatePoint;
                }
            }

            // Try to improve neighbor using owner's nearest wall point
            {
                const Vector candidatePoint = nearestWallPoint_[owner];
                const Scalar newDist =
                    magnitude(neighborCenter - candidatePoint);

                if (newDist < wallDistance_[neighbor])
                {
                    const Scalar change = wallDistance_[neighbor] - newDist;
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
        if (patch.type() != PatchType::wall) continue;

        const BoundaryData& bc =
            bcManager_.fieldBC(patch.patchName(), Field::omega);
        if (bc.type() != BCType::omegaWallFunction) continue;

        wallFunctionFaceIndices_.push_back(faceIdx);

        const size_t cellIdx = face.ownerCell();
        totalWallArea[cellIdx] += face.projectedArea();
        ++wallFaceCountPerCell[cellIdx];
    }

    // Compute per-face weight = faceArea / totalWallArea
    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const size_t cellIdx = face.ownerCell();
        const Scalar area = totalWallArea[cellIdx];

        if (area > S(0.0))
        {
            wallFaceWeight_[face.idx()] = face.projectedArea() / area;
        }

        // cache the owner-cell wall-normal distance
        y_[face.idx()] = std::max
        (
            std::abs(dot(face.dPf(), face.normal())),
            vSmallValue
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
                && patch->get().type() == PatchType::wall)
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
        const size_t cellIdx = wallCellIndices_[i];
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


void kOmegaSST::updateYPlus()
{
    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const size_t cellIdx = face.ownerCell();

        yPlus_[face.idx()] =
            Cmu25_ * std::sqrt(k_[cellIdx]) * y_[face.idx()] / nu_;
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
            const Scalar nutw =
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


void kOmegaSST::strainRate(const TensorField& gradU)
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // S = sqrt(2 * S_ij * S_ij) where S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
        const Scalar symmMagSq = gradU[cellIdx].symm().magnitudeSquared();
        strainRateMag_[cellIdx] = std::sqrt(S(2.0) * symmMagSq);
    }
}


void kOmegaSST::divU(const FaceFluxField& flowRateFace)
{
    const size_t numCells = mesh_.numCells();

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

        divUField_[cellIdx] = sum / cell.volume();
    }
}


void kOmegaSST::kProduction()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar S2 =
            strainRateMag_[cellIdx] * strainRateMag_[cellIdx];

        // k production = nut * S² (unlimited)
        Pk_[cellIdx] = nut_[cellIdx] * S2;
    }
}


void kOmegaSST::updateOmegaWallValues()
{
    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const size_t cellIdx = face.ownerCell();

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
    std::fill(omegaAccum_.begin(), omegaAccum_.end(), S(0.0));

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const size_t cellIdx = face.ownerCell();
        const Scalar faceWeight = wallFaceWeight_[face.idx()];

        if (faceWeight <= S(0.0)) continue;

        if (std::isfinite(omegaWall_[face.idx()]))
        {
            omegaAccum_[cellIdx] += faceWeight * omegaWall_[face.idx()];
        }
    }

    for (size_t i = 0; i < wallCellIndices_.size(); ++i)
    {
        const size_t cellIdx = wallCellIndices_[i];
        wallCellOmega_[i] = std::max(omegaAccum_[cellIdx], smallValue);

        const Scalar f = wallCellFraction_[i];
        omega_[cellIdx] = std::lerp(omega_[cellIdx], wallCellOmega_[i], f);
    }
}


void kOmegaSST::overrideWallCellProduction
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz
)
{
    // Reset per-cell accumulators (capacity from default ctor / initialize())
    wallProductionAccum_.setAll(S(0.0));
    std::fill(hasWallOverride_.begin(), hasWallOverride_.end(), 0);

    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const size_t cellIdx = face.ownerCell();

        if (wallFaceWeight_[face.idx()] <= S(0.0))
        {
            continue;
        }

        if (yPlus_[face.idx()] < yPlusLam_)
        {
            // Viscous sublayer: contribute interior G
            wallProductionAccum_[cellIdx] +=
                wallFaceWeight_[face.idx()] * Pk_[cellIdx];
            hasWallOverride_[cellIdx] = 1;
        }
        else
        {
            // Log layer: G = sqr(uStar*magGradUw*y/uPlus) / (nu*kappa*yPlus)
            const Vector Ucell(Ux[cellIdx], Uy[cellIdx], Uz[cellIdx]);
            const Scalar magGradUw =
                magnitude(Ucell) / y_[face.idx()];
            const Scalar uStar = Cmu25_ * std::sqrt(k_[cellIdx]);
            const Scalar uPlus =
                std::log(std::max(coeffs_.E * yPlus_[face.idx()], S(1.0)))
              / coeffs_.kappa;
            const Scalar uTau2 = uStar * magGradUw * y_[face.idx()] / uPlus;
            const Scalar GWall =
                uTau2 * uTau2
              / (nu_ * coeffs_.kappa * yPlus_[face.idx()]);

            wallProductionAccum_[cellIdx] +=
                wallFaceWeight_[face.idx()] * GWall;
            hasWallOverride_[cellIdx] = 1;
        }
    }

    for (size_t i = 0; i < wallCellIndices_.size(); ++i)
    {
        const size_t cellIdx = wallCellIndices_[i];
        if (!hasWallOverride_[cellIdx]) continue;

        // Blend wall production with interior production using
        // wallCellFraction (wall area / total boundary area)
        const Scalar f = wallCellFraction_[i];
        Pk_[cellIdx] =
            std::lerp(Pk_[cellIdx], wallProductionAccum_[cellIdx], f);
    }
}


void kOmegaSST::crossDiffusion()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        CDkOmega_[cellIdx] =
            S(2.0) * coeffs_.sigmaOmega2
          * dot(gradK_[cellIdx], gradOmega_[cellIdx])
          / std::max(omega_[cellIdx], smallValue);
    }
}


void kOmegaSST::F1()
{
    const size_t numCells = mesh_.numCells();
    constexpr Scalar CDkOmegaMin = S(1e-10);

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar y = std::max(wallDistance_[cellIdx], vSmallValue);
        const Scalar sqrtK = std::sqrt(k_[cellIdx]);

        // Cross-diffusion clipped to positive
        const Scalar CDkw = std::max(CDkOmega_[cellIdx], CDkOmegaMin);

        const Scalar arg1 =
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
        const Scalar arg1Sq = arg1 * arg1;
        f1_[cellIdx] = std::tanh(arg1Sq * arg1Sq);
    }
}


void kOmegaSST::F2()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar y = std::max(wallDistance_[cellIdx], vSmallValue);
        const Scalar sqrtK = std::sqrt(k_[cellIdx]);

        const Scalar arg2 =
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
        f2_[cellIdx] = std::tanh(arg2 * arg2);
    }
}


void kOmegaSST::F3()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar y = std::max(wallDistance_[cellIdx], vSmallValue);

        const Scalar arg3 =
            std::min
            (
                S(150.0) * nu_ / (omega_[cellIdx] * y * y),
                S(10.0)
            );

        // F3 = 1 - tanh(arg3^4)
        const Scalar arg3Sq = arg3 * arg3;
        f3_[cellIdx] = S(1.0) - std::tanh(arg3Sq * arg3Sq);
    }
}


void kOmegaSST::F23()
{
    const size_t numCells = mesh_.numCells();

    if (useF3_)
    {
        #pragma omp parallel for schedule(static)
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            f23_[cellIdx] = f2_[cellIdx] * f3_[cellIdx];
        }
    }
    else
    {
        #pragma omp parallel for schedule(static)
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            f23_[cellIdx] = f2_[cellIdx];
        }
    }
}


void kOmegaSST::omegaProduction()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar S2 =
            strainRateMag_[cellIdx] * strainRateMag_[cellIdx];

        // omega production = gamma * GbyNut (unlimited)
        POmega_[cellIdx] =
            blend(f1_[cellIdx], coeffs_.gamma1, coeffs_.gamma2) * S2;
    }
}


void kOmegaSST::limitProduction()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // Limit k production:
        const Scalar kLimit =
            coeffs_.c1 * coeffs_.betaStar * k_[cellIdx] * omega_[cellIdx];
        Pk_[cellIdx] = std::min(Pk_[cellIdx], kLimit);

        // Limit omega production (Menter 2003 SST)
        const Scalar omegaLimit =
            (coeffs_.c1 / coeffs_.a1) * coeffs_.betaStar * omega_[cellIdx]
          * blend(f1_[cellIdx], coeffs_.gamma1, coeffs_.gamma2)
          * std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23_[cellIdx] * strainRateMag_[cellIdx]
            );
        POmega_[cellIdx] = std::min(POmega_[cellIdx], omegaLimit);
    }
}


void kOmegaSST::solveOmegaEquation(const FaceFluxField& flowRateFace)
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar sigmaOmega =
            blend(f1_[cellIdx], coeffs_.sigmaOmega1, coeffs_.sigmaOmega2);

        GammaOmega_[cellIdx] = nu_ + sigmaOmega * nut_[cellIdx];
    }

    TransportEquation equationOmega
    {
        .field      = Field::omega,
        .phi        = omega_,
        .flowRate   = std::cref(flowRateFace),
        .convScheme = std::cref(omegaConvectionScheme_),
        .Gamma      = std::cref(GammaOmega_),
        .GammaFace  = std::nullopt,
        .source     = omegaSource_,
        .gradPhi    = gradOmega_,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_->buildMatrix(equationOmega);

    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    // Add source terms
    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar cellVolume = mesh_.cells()[cellIdx].volume();

        // Production term: Pk * gamma / nut
        vectorB(eIdx(cellIdx)) += POmega_[cellIdx] * cellVolume;

        // Destruction term: -β·ω² (implicit: β·ω on diagonal)
        const Scalar beta = blend(f1_[cellIdx], coeffs_.beta1, coeffs_.beta2);
        matrixA.coeffRef
        (
            eIdx(cellIdx),
            eIdx(cellIdx)
        ) += beta * omega_[cellIdx] * cellVolume;

        // Cross-diffusion: linearization of (1-F1)*CDkOmega
        const Scalar CDkOmegaLineared =
            (S(1.0) - f1_[cellIdx]) * CDkOmega_[cellIdx]
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
        const Scalar gamma =
            blend(f1_[cellIdx], coeffs_.gamma1, coeffs_.gamma2);
        const Scalar suspOmega =
            (S(2.0) / S(3.0)) * gamma * divUField_[cellIdx];

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

    omegaSolver_->solve(omegaSolution, matrixA, vectorB);

    if (debug_)
    {
        const SolvePerformance& omegaPerformance =
            omegaSolver_->lastPerformance();

        Logger::residualRow
        (
            "omega",
            omegaPerformance.solverName,
            omegaPerformance.iterations,
            omegaPerformance.finalResidual
        );
    }
}


void kOmegaSST::solveKEquation(const FaceFluxField& flowRateFace)
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar sigmaK =
            blend(f1_[cellIdx], coeffs_.sigmaK1, coeffs_.sigmaK2);

        GammaK_[cellIdx] = nu_ + sigmaK * nut_[cellIdx];
    }

    TransportEquation equationK
    {
        .field      = Field::k,
        .phi        = k_,
        .flowRate   = std::cref(flowRateFace),
        .convScheme = std::cref(kConvectionScheme_),
        .Gamma      = std::cref(GammaK_),
        .GammaFace  = std::nullopt,
        .source     = kSource_,
        .gradPhi    = gradK_,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_->buildMatrix(equationK);

    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    // Add k source terms
    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar cellVolume = mesh_.cells()[cellIdx].volume();

        vectorB(eIdx(cellIdx)) += Pk_[cellIdx] * cellVolume;

        // Destruction term: -β*·kω
        const Scalar destruction = coeffs_.betaStar * omega_[cellIdx];
        matrixA.coeffRef(eIdx(cellIdx),eIdx(cellIdx)) +=
            destruction * cellVolume;

        // -(2/3)*divU SuSp term (continuity correction)
        const Scalar suspK = (S(2.0) / S(3.0)) * divUField_[cellIdx];

        matrixA.coeffRef(eIdx(cellIdx),eIdx(cellIdx)) +=
            std::max(suspK, S(0.0)) * cellVolume;

        vectorB(eIdx(cellIdx)) +=
            std::max(-suspK, S(0.0)) * k_[cellIdx] * cellVolume;
    }

    // Apply implicit under-relaxation (Patankar's method)
    matrixConstruct_->relax(alphaK_, k_);

    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    kSolution(k_.data(), eIdx(numCells));

    kSolver_->solve(kSolution, matrixA, vectorB);

    if (debug_)
    {
        const SolvePerformance& kPerformance = kSolver_->lastPerformance();

        Logger::residualRow
        (
            "k",
            kPerformance.solverName,
            kPerformance.iterations,
            kPerformance.finalResidual
        );
    }
}


void kOmegaSST::boundOmega()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar omegaLowerBound =
            k_[cellIdx] / (maxViscosityRatio_ * nu_ + vSmallValue);

        omega_[cellIdx] =
            std::max(omega_[cellIdx], std::max(smallValue, omegaLowerBound));
    }
}


void kOmegaSST::boundK()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        k_[cellIdx] = std::max(k_[cellIdx], smallValue);
    }
}


void kOmegaSST::updateTurbulentViscosity()
{
    // SST turbulent viscosity:
    // nut = a1*k / max(a1*omega, b1*F23*sqrt(S2))
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nut_[cellIdx] =
            (coeffs_.a1 * k_[cellIdx])
          / std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23_[cellIdx] * strainRateMag_[cellIdx]
            );
    }
}


void kOmegaSST::updateWallShearStress
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz
)
{
    for (size_t faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const size_t cellIdx = face.ownerCell();

        // Project velocity onto wall plane (tangential component)
        const Vector Ucell(Ux[cellIdx], Uy[cellIdx], Uz[cellIdx]);
        const Scalar normalVelocity = dot(Ucell, face.normal());
        const Vector tangentVelocity =
            Ucell - face.normal() * normalVelocity;
        const Scalar tangentVelocityMag = magnitude(tangentVelocity);

        const Scalar tau =
            yPlus_[face.idx()] < yPlusLam_
          ? nu_ * tangentVelocityMag / y_[face.idx()]
          : Cmu25_ * Cmu25_ * k_[cellIdx];

        wallShearStress_[face.idx()] = std::min(tau, S(1000.0));
    }
}


void kOmegaSST::logFieldDiagnostics() const
{
    if (!debug_) return;

    const size_t numCells = mesh_.numCells();
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

    const Scalar n = S(numCells);

    Logger::subsection("Turbulence field statistics");
    Logger::scalarStat("k",     kMin, kMax, kSum / n);
    Logger::scalarStat("omega", oMin, oMax, oSum / n);
    Logger::scalarStat("nut",   nMin, nMax, nSum / n);
}
