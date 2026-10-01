/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file kOmegaSST.cpp
 * @brief Implementation of k-omega SST turbulence model
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "kOmegaSST.h"

// Standard library headers
#include <cmath>
#include <algorithm>
#include <limits>

// Project headers
#include "HaloExchange.h"
#include "Logger.h"
#include "Matrix.h"
#include "Reduce.h"
#include "BoundaryConditions.h"
#include "GradientScheme.h"
#include "ConvectionScheme.h"
#include "LinearSolvers.h"

// ************************* Special Member Functions *************************

kOmegaSST::kOmegaSST
(
    const Mesh& mesh,
    const BoundaryConditions& bc,
    const TimeScheme& timeScheme,
    const GradientScheme& gradientScheme,
    const ConvectionScheme& kScheme,
    LinearSolver& kSolver,
    const ConvectionScheme& omegaScheme,
    LinearSolver& omegaSolver,
    Scalar deltaT,
    Scalar nu,
    Scalar initialK,
    Scalar initialOmega,
    Scalar alphaK,
    Scalar alphaOmega,
    bool roughWall,
    bool debug
)
:
    RANS
    {
        mesh,
        bc,
        timeScheme,
        gradientScheme,
        kScheme,
        kSolver,
        omegaScheme,
        omegaSolver,
        deltaT,
        nu,
        alphaK,
        alphaOmega,
        debug
    }
{
    useF3_ = roughWall;
    // Compute yPlusLam and wall-function geometry
    updateYPlusLam(coeffs_.kappa, coeffs_.E);
    updateWallDistance();
    initializeWallFunctionGeometry(bcManager(), Field::omega);
    wallCellOmega_.assign(wallCellIndices().size(), S(0.0));

    // The ghosts show the neighbor's constrained cells
    wallConstraintFraction_.setAll(S(0.0));

    for (Index i = 0; i < wallCellIndices().size(); ++i)
    {
        wallConstraintFraction_[wallCellIndices()[i]] =
            wallCellFraction()[i];
    }

    Halo::exchange({&wallConstraintFraction_});

    // Initialize turbulence fields with initial conditions
    k().setAll(initialK);
    omega_.setAll(initialOmega);

    // nut = k/omega until strain rate and F23 reach the SST limiter
    const Count numCells = mesh.numCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nut()[cellIdx] =
            std::max(k()[cellIdx] / (omega_[cellIdx] + vSmallValue), S(0.0));
    }

    updateYPlus();
    updateNutWall();
}

kOmegaSST::~kOmegaSST() noexcept = default;

// ****************************** Solve kOmegaSST *****************************

void kOmegaSST::solve
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz,
    const FaceFluxField& flowRateFace,
    const TensorField& gradU
)
{
    // Snapshot fields before the solve for residual computation
    kPrev() = k();
    omegaPrev_ = omega_;

    // Update y+ on wall faces
    updateYPlus();

    // Compute geometric quantities
    const ScalarField strainRateSq = strainRateSquared(gradU);
    const ScalarField divU = velocityDivergence(flowRateFace);

    // Compute k Production
    ScalarField Pk = kProduction(strainRateSq);

    // Update wall-function boundary values for omega
    updateOmegaWallValues();

    // Pre-set wall-cell omega via area-weighted lerp
    applyOmegaWallCellValues();
    Halo::exchange({&omega_});

    // Override k production at wall-adjacent cells
    overrideWallCellProduction(Ux, Uy, Uz, Pk);

    // Compute gradients and cross-diffusion
    gradientScheme().fieldGradient(Field::k, k(), gradK());
    VectorField gradOmega(mesh());
    gradientScheme().fieldGradient(Field::omega, omega_, gradOmega);

    // The k/omega deferred corrections read both gradients at cut faces
    Halo::exchange({&gradK(), &gradOmega});

    const ScalarField CDkOmega = crossDiffusion(gradOmega);

    // Compute blending functions
    const ScalarField f1 = blendingF1(CDkOmega);
    const ScalarField f2 = blendingF2();
    const ScalarField f23 =
        useF3_ ? blendingF23(f2, blendingF3()) : f2;

    // Compute omega production
    ScalarField POmega = omegaProduction(f1, strainRateSq);

    // Apply SST production limiters
    limitProduction(f1, f23, strainRateSq, Pk, POmega);

    // Solve omega transport equation
    solveOmegaEquation(flowRateFace, divU, f1, CDkOmega, POmega, gradOmega);
    boundOmega();
    Halo::exchange({&omega_});

    // Solve k transport equation
    solveKEquation(flowRateFace, divU, f1, Pk);
    boundK();
    boundOmega();
    Halo::exchange({&k(), &omega_});

    // Update turbulent viscosity with SST limiter
    nut() = computeTurbulentViscosity(f23, strainRateSq);

    // Update wall-function nut on wall faces
    updateNutWall();
    Halo::exchange({&nut()});

    // Compute normalised k/omega change against the pre-solve snapshots
    updateResiduals(omega_, omegaPrev_);

    // Log min/max/mean of k, omega, nut
    const Count numCells = mesh().numDomainCells();

    if (debug())
    {
        Scalar kMin = std::numeric_limits<Scalar>::max();
        Scalar kMax = std::numeric_limits<Scalar>::lowest();
        Scalar kSum = S(0.0);

        Scalar omegaMin = std::numeric_limits<Scalar>::max();
        Scalar omegaMax = std::numeric_limits<Scalar>::lowest();
        Scalar omegaSum = S(0.0);

        Scalar nutMin = std::numeric_limits<Scalar>::max();
        Scalar nutMax = std::numeric_limits<Scalar>::lowest();
        Scalar nutSum = S(0.0);

        for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            kMin = std::min(kMin, k()[cellIdx]);
            kMax = std::max(kMax, k()[cellIdx]);
            kSum += k()[cellIdx];

            omegaMin = std::min(omegaMin, omega_[cellIdx]);
            omegaMax = std::max(omegaMax, omega_[cellIdx]);
            omegaSum += omega_[cellIdx];

            nutMin = std::min(nutMin, nut()[cellIdx]);
            nutMax = std::max(nutMax, nut()[cellIdx]);
            nutSum += nut()[cellIdx];
        }

        kMin = globalMin(kMin);
        kMax = globalMax(kMax);
        kSum = globalSum(kSum);

        omegaMin = globalMin(omegaMin);
        omegaMax = globalMax(omegaMax);
        omegaSum = globalSum(omegaSum);

        nutMin = globalMin(nutMin);
        nutMax = globalMax(nutMax);
        nutSum = globalSum(nutSum);

        const Scalar n = S(globalSum(numCells));

        Logger::subsection("Turbulence field statistics");
        Logger::scalarStat("k", kMin, kMax, kSum / n);
        Logger::scalarStat("omega", omegaMin, omegaMax, omegaSum / n);
        Logger::scalarStat("nut", nutMin, nutMax, nutSum / n);
    }
}

// ************************ Inlet Condition Calculators ***********************

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

void kOmegaSST::updateNutWall()
{
    nutWall().setAll(S(0.0));

    for (Index i = 0; i < wallFunctionFaceIndices().size(); ++i)
    {
        const Index faceIdx = wallFunctionFaceIndices()[i];
        const auto& face = mesh().faces()[faceIdx];

        if (yPlus()[face.idx()] > yPlusLam())
        {
            // Log layer: nutw = nu * (yPlus*kappa/ln(E*yPlus) - 1)
            const Scalar nutw =
                nu()
              * (
                    yPlus()[face.idx()] * coeffs_.kappa
                  / std::log(std::max(coeffs_.E * yPlus()[face.idx()], S(1.0)))
                  - S(1.0)
                );

            nutWall()[face.idx()] = std::max(nutw, S(0.0));
        }
        // Viscous sublayer: nutWall = 0 (already initialized)
    }
}


void kOmegaSST::updateOmegaWallValues()
{
    for (Index i = 0; i < wallFunctionFaceIndices().size(); ++i)
    {
        const Index faceIdx = wallFunctionFaceIndices()[i];
        const auto& face = mesh().faces()[faceIdx];
        const Index cellIdx = face.ownerCell();

        if (yPlus()[face.idx()] < yPlusLam())
        {
            omegaWall_[face.idx()] =
                S(6.0) * nu()
              / (coeffs_.beta1 * y()[face.idx()] * y()[face.idx()]);
        }
        else
        {
            omegaWall_[face.idx()] =
                std::sqrt(k()[cellIdx])
              / (Cmu25_ * coeffs_.kappa * y()[face.idx()]);
        }
    }
}


void kOmegaSST::applyOmegaWallCellValues()
{
    ScalarList omegaAccum(mesh().numCells(), S(0.0));

    for (Index faceIdx : wallFunctionFaceIndices())
    {
        const auto& face = mesh().faces()[faceIdx];
        const Index cellIdx = face.ownerCell();
        const Scalar faceWeight = wallFaceWeight()[face.idx()];

        if (faceWeight <= S(0.0)) continue;

        if (std::isfinite(omegaWall_[face.idx()]))
        {
            omegaAccum[cellIdx] += faceWeight * omegaWall_[face.idx()];
        }
    }

    for (Index i = 0; i < wallCellIndices().size(); ++i)
    {
        const Index cellIdx = wallCellIndices()[i];
        wallCellOmega_[i] = std::max(omegaAccum[cellIdx], smallValue);

        const Scalar f = wallCellFraction()[i];
        omega_[cellIdx] = std::lerp(omega_[cellIdx], wallCellOmega_[i], f);
    }
}


void kOmegaSST::overrideWallCellProduction
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz,
    ScalarField& Pk
)
{
    ScalarField wallProductionAccum(mesh(), S(0.0));
    std::vector<char> hasWallOverride(mesh().numCells(), 0);

    for (Index faceIdx : wallFunctionFaceIndices())
    {
        const auto& face = mesh().faces()[faceIdx];
        const Index cellIdx = face.ownerCell();

        if (wallFaceWeight()[face.idx()] <= S(0.0))
        {
            continue;
        }

        if (yPlus()[face.idx()] < yPlusLam())
        {
            // Viscous sublayer: contribute interior G
            wallProductionAccum[cellIdx] +=
                wallFaceWeight()[face.idx()] * Pk[cellIdx];
            hasWallOverride[cellIdx] = 1;
        }
        else
        {
            // Log layer: G = sqr(uStar*magGradUw*y/uPlus) / (nu*kappa*yPlus)
            const Vector Ucell(Ux[cellIdx], Uy[cellIdx], Uz[cellIdx]);
            const Scalar magGradUw =
                magnitude(Ucell) / y()[face.idx()];
            const Scalar uStar = Cmu25_ * std::sqrt(k()[cellIdx]);
            const Scalar uPlus =
                std::log(std::max(coeffs_.E * yPlus()[face.idx()], S(1.0)))
              / coeffs_.kappa;
            const Scalar uTau2 = uStar * magGradUw * y()[face.idx()] / uPlus;
            const Scalar GWall =
                uTau2 * uTau2
              / (nu() * coeffs_.kappa * yPlus()[face.idx()]);

            wallProductionAccum[cellIdx] +=
                wallFaceWeight()[face.idx()] * GWall;
            hasWallOverride[cellIdx] = 1;
        }
    }

    for (Index i = 0; i < wallCellIndices().size(); ++i)
    {
        const Index cellIdx = wallCellIndices()[i];
        if (!hasWallOverride[cellIdx]) continue;

        // Blend by wallCellFraction (wall area / total boundary area)
        const Scalar f = wallCellFraction()[i];
        Pk[cellIdx] =
            std::lerp(Pk[cellIdx], wallProductionAccum[cellIdx], f);
    }
}


ScalarField kOmegaSST::kProduction
(
    const ScalarField& strainRateSq
) const
{
    const Count numCells = mesh().numDomainCells();
    ScalarField Pk{mesh()};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // k production = nut * S² (unlimited)
        Pk[cellIdx] = nut()[cellIdx] * strainRateSq[cellIdx];
    }

    return Pk;
}


ScalarField kOmegaSST::crossDiffusion
(
    const VectorField& gradOmega
) const
{
    const Count numCells = mesh().numDomainCells();
    ScalarField CDkOmega{mesh()};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        CDkOmega[cellIdx] =
            S(2.0) * coeffs_.sigmaOmega2
          * dot(gradK()[cellIdx], gradOmega[cellIdx])
          / std::max(omega_[cellIdx], smallValue);
    }

    return CDkOmega;
}


ScalarField kOmegaSST::blendingF1
(
    const ScalarField& CDkOmega
) const
{
    const Count numCells = mesh().numDomainCells();
    constexpr Scalar CDkOmegaMin = S(1e-10);
    ScalarField f1{mesh()};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar y = std::max(wallDistance()[cellIdx], vSmallValue);
        const Scalar sqrtK = std::sqrt(k()[cellIdx]);

        // Cross-diffusion clipped to positive
        const Scalar CDkw = std::max(CDkOmega[cellIdx], CDkOmegaMin);

        const Scalar arg1 =
            std::min
            (
                std::min
                (
                    std::max
                    (
                        sqrtK / (coeffs_.betaStar * omega_[cellIdx] * y),
                        S(500.0) * nu() / (omega_[cellIdx] * y * y)
                    ),
                    S(4.0) * coeffs_.sigmaOmega2 * k()[cellIdx] / (CDkw * y * y)
                ),
                S(10.0)
            );

        // F1 = tanh(arg1^4)
        // Note: Direct multiplication is faster than std::pow(arg1, 4)
        const Scalar arg1Sq = arg1 * arg1;
        f1[cellIdx] = std::tanh(arg1Sq * arg1Sq);
    }

    return f1;
}


ScalarField kOmegaSST::blendingF2() const
{
    const Count numCells = mesh().numDomainCells();
    ScalarField f2{mesh()};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar y = std::max(wallDistance()[cellIdx], vSmallValue);
        const Scalar sqrtK = std::sqrt(k()[cellIdx]);

        const Scalar arg2 =
            std::min
            (
                std::max
                (
                    S(2.0) * sqrtK
                  / (coeffs_.betaStar * omega_[cellIdx] * y),
                    S(500.0) * nu() / (omega_[cellIdx] * y * y)
                ),
                S(100.0)
            );

        // F2 = tanh(arg2^2)
        f2[cellIdx] = std::tanh(arg2 * arg2);
    }

    return f2;
}


ScalarField kOmegaSST::blendingF3() const
{
    const Count numCells = mesh().numDomainCells();
    ScalarField f3{mesh()};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar y = std::max(wallDistance()[cellIdx], vSmallValue);

        const Scalar arg3 =
            std::min
            (
                S(150.0) * nu() / (omega_[cellIdx] * y * y),
                S(10.0)
            );

        // F3 = 1 - tanh(arg3^4)
        const Scalar arg3Sq = arg3 * arg3;
        f3[cellIdx] = S(1.0) - std::tanh(arg3Sq * arg3Sq);
    }

    return f3;
}


ScalarField kOmegaSST::blendingF23
(
    const ScalarField& f2,
    const ScalarField& f3
) const
{
    const Count numCells = mesh().numDomainCells();
    ScalarField f23{mesh()};

    if (useF3_)
    {
        for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            f23[cellIdx] = f2[cellIdx] * f3[cellIdx];
        }
    }
    else
    {
        for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            f23[cellIdx] = f2[cellIdx];
        }
    }

    return f23;
}


ScalarField kOmegaSST::omegaProduction
(
    const ScalarField& f1,
    const ScalarField& strainRateSq
) const
{
    const Count numCells = mesh().numDomainCells();
    ScalarField POmega{mesh()};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // omega production = gamma * S² (unlimited)
        POmega[cellIdx] =
            blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2)
          * strainRateSq[cellIdx];
    }

    return POmega;
}


ScalarField kOmegaSST::Gamma
(
    const ScalarField& f1,
    Scalar sigma1,
    Scalar sigma2
) const
{
    const Count numCells = mesh().numDomainCells();
    ScalarField Gamma{mesh()};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar sigma = blend(f1[cellIdx], sigma1, sigma2);
        Gamma[cellIdx] = nu() + sigma * nut()[cellIdx];
    }

    return Gamma;
}


void kOmegaSST::limitProduction
(
    const ScalarField& f1,
    const ScalarField& f23,
    const ScalarField& strainRateSq,
    ScalarField& Pk,
    ScalarField& POmega
) const
{
    const Count numCells = mesh().numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // Limit k production:
        const Scalar kLimit =
            coeffs_.c1 * coeffs_.betaStar * k()[cellIdx] * omega_[cellIdx];
        Pk[cellIdx] = std::min(Pk[cellIdx], kLimit);

        // Limit omega production (Menter 2003 SST)
        const Scalar strainRate = std::sqrt(strainRateSq[cellIdx]);
        const Scalar omegaLimit =
            (coeffs_.c1 / coeffs_.a1) * coeffs_.betaStar * omega_[cellIdx]
          * blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2)
          * std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23[cellIdx] * strainRate
            );
        POmega[cellIdx] = std::min(POmega[cellIdx], omegaLimit);
    }
}


void kOmegaSST::solveOmegaEquation
(
    const FaceFluxField& flowRateFace,
    const ScalarField& divU,
    const ScalarField& f1,
    const ScalarField& CDkOmega,
    const ScalarField& POmega,
    const VectorField& gradOmega
)
{
    const Count numCells = mesh().numDomainCells();
    ScalarField GammaOmega =
        Gamma(f1, coeffs_.sigmaOmega1, coeffs_.sigmaOmega2);
    Halo::exchange({&GammaOmega});
    cellToFaceDiffusion(GammaOmega, gammaOmegaFace_);

    const ScalarField omegaSource{mesh(), S(0.0)};

    TransportEquation equationOmega
    {
        .field      = Field::omega,
        .phi        = omega_,
        .transient  =
            ddtTerm(dissipationPrevStep(), dissipationDdtPrevStep()),
        .convection =
            ConvectionTerm{flowRateFace, dissipationConvectionScheme()},
        .GammaFace  = gammaOmegaFace_,
        .source     = omegaSource,
        .gradPhi    = gradOmega,
        .gradScheme = gradientScheme()
    };

    matrixConstruct()->buildMatrix(equationOmega);

    const std::span<Scalar> diagonal = matrixConstruct()->diagonal();
    ScalarList& vectorB = matrixConstruct()->vectorB();

    // Add source terms
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar cellVolume = mesh().cells()[cellIdx].volume();

        // Production term: add the limited omega production POmega to RHS
        vectorB[cellIdx] += POmega[cellIdx] * cellVolume;

        // Destruction term: -β·ω² (implicit: β·ω on diagonal)
        const Scalar beta = blend(f1[cellIdx], coeffs_.beta1, coeffs_.beta2);
        diagonal[cellIdx] += beta * omega_[cellIdx] * cellVolume;

        // Cross-diffusion: linearization of (1-F1)*CDkOmega
        const Scalar CDkOmegaLineared =
            (S(1.0) - f1[cellIdx]) * CDkOmega[cellIdx]
          / (omega_[cellIdx] + vSmallValue);

        if (CDkOmegaLineared < S(0.0))
        {
            diagonal[cellIdx] += -CDkOmegaLineared * cellVolume;
        }
        else
        {
            vectorB[cellIdx] +=
                CDkOmegaLineared * omega_[cellIdx] * cellVolume;
        }

        // -(2/3)*gamma*divU SuSp term (continuity correction)
        const Scalar gamma =
            blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2);
        const Scalar suspOmega =
            (S(2.0) / S(3.0)) * gamma * divU[cellIdx];

        diagonal[cellIdx] += std::max(suspOmega, S(0.0)) * cellVolume;

        vectorB[cellIdx] +=
            std::max(-suspOmega, S(0.0)) * omega_[cellIdx] * cellVolume;
    }

    // Apply under-relaxation
    matrixConstruct()->relax(alphaDissipation(), omega_);

    // Fix wall-cell rows to impose omega = omegaWall
    matrixConstruct()->setValues
    (
        wallCellIndices(),
        wallCellOmega_,
        wallConstraintFraction_,
        omega_,
        wallCellFraction()
    );

    matrixConstruct()->assemble();

    dissipationSolver().solve
    (
        {omega_.data(), numCells},
        matrixConstruct()->matrixA(),
        matrixConstruct()->rhsVec()
    );

    if (debug())
    {
        const SolvePerformance& omegaPerformance =
            dissipationSolver().lastPerformance();

        Logger::residualRow
        (
            "omega",
            omegaPerformance.solverName,
            omegaPerformance.iterations,
            omegaPerformance.finalResidual
        );
    }
}


void kOmegaSST::boundOmega()
{
    const Count numCells = mesh().numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar omegaLowerBound =
            k()[cellIdx] / (maxViscosityRatio_ * nu() + vSmallValue);

        omega_[cellIdx] =
            std::max(omega_[cellIdx], std::max(smallValue, omegaLowerBound));
    }
}


void kOmegaSST::solveKEquation
(
    const FaceFluxField& flowRateFace,
    const ScalarField& divU,
    const ScalarField& f1,
    const ScalarField& Pk
)
{
    const Count numCells = mesh().numDomainCells();
    ScalarField GammaK =
        Gamma(f1, coeffs_.sigmaK1, coeffs_.sigmaK2);
    Halo::exchange({&GammaK});
    cellToFaceDiffusion(GammaK, gammaKFace_);

    const ScalarField kSource{mesh(), S(0.0)};

    TransportEquation equationK
    {
        .field      = Field::k,
        .phi        = k(),
        .transient  = ddtTerm(kPrevStep(), kDdtPrevStep()),
        .convection = ConvectionTerm{flowRateFace, kConvectionScheme()},
        .GammaFace  = gammaKFace_,
        .source     = kSource,
        .gradPhi    = gradK(),
        .gradScheme = gradientScheme()
    };

    matrixConstruct()->buildMatrix(equationK);

    const std::span<Scalar> diagonal = matrixConstruct()->diagonal();
    ScalarList& vectorB = matrixConstruct()->vectorB();

    // Add k source terms
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar cellVolume = mesh().cells()[cellIdx].volume();

        vectorB[cellIdx] += Pk[cellIdx] * cellVolume;

        // Destruction term: -β*·kω
        const Scalar destruction = coeffs_.betaStar * omega_[cellIdx];
        diagonal[cellIdx] += destruction * cellVolume;

        // -(2/3)*divU SuSp term (continuity correction)
        const Scalar suspK = (S(2.0) / S(3.0)) * divU[cellIdx];

        diagonal[cellIdx] += std::max(suspK, S(0.0)) * cellVolume;

        vectorB[cellIdx] +=
            std::max(-suspK, S(0.0)) * k()[cellIdx] * cellVolume;
    }

    // Apply implicit under-relaxation (Patankar's method)
    matrixConstruct()->relax(alphaK(), k());

    matrixConstruct()->assemble();

    kSolver().solve
    (
        {k().data(), numCells},
        matrixConstruct()->matrixA(),
        matrixConstruct()->rhsVec()
    );

    if (debug())
    {
        const SolvePerformance& kPerformance = kSolver().lastPerformance();

        Logger::residualRow
        (
            "k",
            kPerformance.solverName,
            kPerformance.iterations,
            kPerformance.finalResidual
        );
    }
}


void kOmegaSST::boundK()
{
    const Count numCells = mesh().numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        k()[cellIdx] = std::max(k()[cellIdx], smallValue);
    }
}


ScalarField kOmegaSST::computeTurbulentViscosity
(
    const ScalarField& f23,
    const ScalarField& strainRateSq
) const
{
    // SST turbulent viscosity:
    // nut = a1*k / max(a1*omega, F23*sqrt(S2))
    const Count numCells = mesh().numDomainCells();
    ScalarField nut{mesh()};

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar strainRate = std::sqrt(strainRateSq[cellIdx]);
        nut[cellIdx] =
            (coeffs_.a1 * k()[cellIdx])
          / std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23[cellIdx] * strainRate
            );
    }

    return nut;
}
