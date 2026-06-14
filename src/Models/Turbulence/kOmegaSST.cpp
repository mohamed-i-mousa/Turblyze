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

// External library headers
#include <omp.h>

// Project headers
#include "Logger.h"
#include "Matrix.h"
#include "BoundaryConditions.h"
#include "GradientScheme.h"
#include "ConvectionSchemes.h"
#include "LinearSolvers.h"

// ************************* Special Member Functions *************************

kOmegaSST::kOmegaSST
(
    const Mesh& mesh,
    const BoundaryConditions& bc,
    const GradientScheme& gradientScheme,
    const ConvectionSchemes& kScheme,
    LinearSolver& kSolver,
    const ConvectionSchemes& omegaScheme,
    LinearSolver& omegaSolver,
    Scalar nu,
    Scalar initialK,
    Scalar initialOmega,
    Scalar alphaK,
    Scalar alphaOmega,
    bool debug
)
:
    RANS
    {
        mesh,
        bc,
        gradientScheme,
        kScheme,
        kSolver,
        omegaScheme,
        omegaSolver,
        nu,
        alphaK,
        alphaOmega,
        debug
    }
{
    // Compute yPlusLam and wall-function geometry
    updateYPlusLam(coeffs_.kappa, coeffs_.E);
    updateWallDistance();
    initializeWallFunctionGeometry
    (
        bcManager_,
        Field::omega,
        BCType::omegaWallFunction
    );
    wallCellOmega_.assign(wallCellIndices_.size(), S(0.0));

    // Initialize turbulence fields with initial conditions
    k_.setAll(initialK);
    omega_.setAll(initialOmega);

    // initialize nut = k/omega
    initializeTurbulentViscosity();
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
    kPrev_ = k_;
    omegaPrev_ = omega_;

    // Update y+ on wall faces
    updateYPlus();

    // Compute geometric quantities
    const ScalarField strainRateMag = computeStrainRateMagnitude(gradU);
    const ScalarField divU = velocityDivergence(flowRateFace);

    // Compute k Production
    ScalarField Pk = kProduction(strainRateMag);

    // Update wall-function boundary values for omega
    updateOmegaWallValues();

    // Pre-set wall-cell omega via area-weighted lerp
    applyOmegaWallCellValues();

    // Override k production at wall-adjacent cells
    overrideWallCellProduction(Ux, Uy, Uz, Pk);

    // Compute gradients and cross-diffusion
    gradientScheme_.fieldGradient(Field::k, k_, gradK_);
    VectorField gradOmega;
    gradientScheme_.fieldGradient(Field::omega, omega_, gradOmega);

    const ScalarField CDkOmega = crossDiffusion(gradOmega);

    // Compute blending functions
    const ScalarField f1 = blendingF1(CDkOmega);
    const ScalarField f2 = blendingF2();
    const ScalarField f3 = blendingF3();
    const ScalarField f23 = blendingF23(f2, f3);

    // Compute omega production
    ScalarField POmega = omegaProduction(f1, strainRateMag);

    // Apply SST production limiters
    limitProduction(f1, f23, strainRateMag, Pk, POmega);

    // Solve omega transport equation
    solveOmegaEquation(flowRateFace, divU, f1, CDkOmega, POmega, gradOmega);
    boundOmega();

    // Solve k transport equation
    solveKEquation(flowRateFace, divU, f1, Pk);
    boundK();
    boundOmega();

    // Update turbulent viscosity with SST limiter
    nut_ = computeTurbulentViscosity(f23, strainRateMag);

    // Update wall-function nut on wall faces
    updateNutWall();

    // Compute normalised k/omega change against the pre-solve snapshots
    updateResiduals(omega_, omegaPrev_);

    // Log min/max/mean of k, omega, nut
    logFieldDiagnostics();
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

void kOmegaSST::initializeTurbulentViscosity()
{
    // Simple k-omega estimate: nut = k / omega
    // Used before strain rate & F23 are available for the SST limiter.
    const Count numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nut_[cellIdx] =
            std::max(k_[cellIdx] / (omega_[cellIdx] + vSmallValue), S(0.0));
    }
}


void kOmegaSST::updateNutWall()
{
    nutWall_.setAll(S(0.0));

    #pragma omp parallel for schedule(static)
    for (Index i = 0; i < wallFunctionFaceIndices_.size(); ++i)
    {
        const Index faceIdx = wallFunctionFaceIndices_[i];
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


void kOmegaSST::updateOmegaWallValues()
{
    #pragma omp parallel for schedule(static)
    for (Index i = 0; i < wallFunctionFaceIndices_.size(); ++i)
    {
        const Index faceIdx = wallFunctionFaceIndices_[i];
        const auto& face = mesh_.faces()[faceIdx];
        const Index cellIdx = face.ownerCell();

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
    ScalarList omegaAccum(mesh_.numCells(), S(0.0));

    for (Index faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const Index cellIdx = face.ownerCell();
        const Scalar faceWeight = wallFaceWeight_[face.idx()];

        if (faceWeight <= S(0.0)) continue;

        if (std::isfinite(omegaWall_[face.idx()]))
        {
            omegaAccum[cellIdx] += faceWeight * omegaWall_[face.idx()];
        }
    }

    for (Index i = 0; i < wallCellIndices_.size(); ++i)
    {
        const Index cellIdx = wallCellIndices_[i];
        wallCellOmega_[i] = std::max(omegaAccum[cellIdx], smallValue);

        const Scalar f = wallCellFraction_[i];
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
    ScalarField wallProductionAccum;
    wallProductionAccum.setAll(S(0.0));
    std::vector<char> hasWallOverride(mesh_.numCells(), 0);

    for (Index faceIdx : wallFunctionFaceIndices_)
    {
        const auto& face = mesh_.faces()[faceIdx];
        const Index cellIdx = face.ownerCell();

        if (wallFaceWeight_[face.idx()] <= S(0.0))
        {
            continue;
        }

        if (yPlus_[face.idx()] < yPlusLam_)
        {
            // Viscous sublayer: contribute interior G
            wallProductionAccum[cellIdx] +=
                wallFaceWeight_[face.idx()] * Pk[cellIdx];
            hasWallOverride[cellIdx] = 1;
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

            wallProductionAccum[cellIdx] +=
                wallFaceWeight_[face.idx()] * GWall;
            hasWallOverride[cellIdx] = 1;
        }
    }

    for (Index i = 0; i < wallCellIndices_.size(); ++i)
    {
        const Index cellIdx = wallCellIndices_[i];
        if (!hasWallOverride[cellIdx]) continue;

        // Blend wall production with interior production using
        // wallCellFraction (wall area / total boundary area)
        const Scalar f = wallCellFraction_[i];
        Pk[cellIdx] =
            std::lerp(Pk[cellIdx], wallProductionAccum[cellIdx], f);
    }
}


ScalarField kOmegaSST::kProduction
(
    const ScalarField& strainRateMag
) const
{
    const Count numCells = mesh_.numCells();
    ScalarField Pk;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar S2 =
            strainRateMag[cellIdx] * strainRateMag[cellIdx];

        // k production = nut * S² (unlimited)
        Pk[cellIdx] = nut_[cellIdx] * S2;
    }

    return Pk;
}


ScalarField kOmegaSST::crossDiffusion
(
    const VectorField& gradOmega
) const
{
    const Count numCells = mesh_.numCells();
    ScalarField CDkOmega;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        CDkOmega[cellIdx] =
            S(2.0) * coeffs_.sigmaOmega2
          * dot(gradK_[cellIdx], gradOmega[cellIdx])
          / std::max(omega_[cellIdx], smallValue);
    }

    return CDkOmega;
}


ScalarField kOmegaSST::blendingF1
(
    const ScalarField& CDkOmega
) const
{
    const Count numCells = mesh_.numCells();
    constexpr Scalar CDkOmegaMin = S(1e-10);
    ScalarField f1;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar y = std::max(wallDistance_[cellIdx], vSmallValue);
        const Scalar sqrtK = std::sqrt(k_[cellIdx]);

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
                        S(500.0) * nu_ / (omega_[cellIdx] * y * y)
                    ),
                    S(4.0) * coeffs_.sigmaOmega2 * k_[cellIdx] / (CDkw * y * y)
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
    const Count numCells = mesh_.numCells();
    ScalarField f2;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
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
        f2[cellIdx] = std::tanh(arg2 * arg2);
    }

    return f2;
}


ScalarField kOmegaSST::blendingF3() const
{
    const Count numCells = mesh_.numCells();
    ScalarField f3;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
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
    const Count numCells = mesh_.numCells();
    ScalarField f23;

    if (useF3_)
    {
        #pragma omp parallel for schedule(static)
        for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            f23[cellIdx] = f2[cellIdx] * f3[cellIdx];
        }
    }
    else
    {
        #pragma omp parallel for schedule(static)
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
    const ScalarField& strainRateMag
) const
{
    const Count numCells = mesh_.numCells();
    ScalarField POmega;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar S2 =
            strainRateMag[cellIdx] * strainRateMag[cellIdx];

        // omega production = gamma * GbyNut (unlimited)
        POmega[cellIdx] =
            blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2) * S2;
    }

    return POmega;
}


ScalarField kOmegaSST::computeGammaK(const ScalarField& f1) const
{
    const Count numCells = mesh_.numCells();
    ScalarField GammaK;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar sigmaK =
            blend(f1[cellIdx], coeffs_.sigmaK1, coeffs_.sigmaK2);

        GammaK[cellIdx] = nu_ + sigmaK * nut_[cellIdx];
    }

    return GammaK;
}


ScalarField kOmegaSST::computeGammaOmega(const ScalarField& f1) const
{
    const Count numCells = mesh_.numCells();
    ScalarField GammaOmega;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar sigmaOmega =
            blend(f1[cellIdx], coeffs_.sigmaOmega1, coeffs_.sigmaOmega2);

        GammaOmega[cellIdx] = nu_ + sigmaOmega * nut_[cellIdx];
    }

    return GammaOmega;
}


void kOmegaSST::limitProduction
(
    const ScalarField& f1,
    const ScalarField& f23,
    const ScalarField& strainRateMag,
    ScalarField& Pk,
    ScalarField& POmega
) const
{
    const Count numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        // Limit k production:
        const Scalar kLimit =
            coeffs_.c1 * coeffs_.betaStar * k_[cellIdx] * omega_[cellIdx];
        Pk[cellIdx] = std::min(Pk[cellIdx], kLimit);

        // Limit omega production (Menter 2003 SST)
        const Scalar omegaLimit =
            (coeffs_.c1 / coeffs_.a1) * coeffs_.betaStar * omega_[cellIdx]
          * blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2)
          * std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23[cellIdx] * strainRateMag[cellIdx]
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
    const Count numCells = mesh_.numCells();
    const ScalarField GammaOmega = computeGammaOmega(f1);
    cellToFaceDiffusion(GammaOmega, gammaOmegaFace_);

    const ScalarField omegaSource{S(0.0)};

    TransportEquation equationOmega
    {
        .field      = Field::omega,
        .phi        = omega_,
        .convection =
            ConvectionTerm{flowRateFace, dissipationConvectionScheme_},
        .GammaFace  = gammaOmegaFace_,
        .source     = omegaSource,
        .gradPhi    = gradOmega,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_->buildMatrix(equationOmega);

    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    // Add source terms
    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar cellVolume = mesh_.cells()[cellIdx].volume();

        // Production term: add the limited omega production POmega to RHS
        vectorB(eIdx(cellIdx)) += POmega[cellIdx] * cellVolume;

        // Destruction term: -β·ω² (implicit: β·ω on diagonal)
        const Scalar beta = blend(f1[cellIdx], coeffs_.beta1, coeffs_.beta2);
        matrixA.coeffRef
        (
            eIdx(cellIdx),
            eIdx(cellIdx)
        ) += beta * omega_[cellIdx] * cellVolume;

        // Cross-diffusion: linearization of (1-F1)*CDkOmega
        const Scalar CDkOmegaLineared =
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
        const Scalar gamma =
            blend(f1[cellIdx], coeffs_.gamma1, coeffs_.gamma2);
        const Scalar suspOmega =
            (S(2.0) / S(3.0)) * gamma * divU[cellIdx];

        matrixA.coeffRef
        (
            eIdx(cellIdx),
            eIdx(cellIdx)
        ) += std::max(suspOmega, S(0.0)) * cellVolume;

        vectorB(eIdx(cellIdx)) +=
            std::max(-suspOmega, S(0.0)) * omega_[cellIdx] * cellVolume;
    }

    // Apply under-relaxation
    matrixConstruct_->relax(alphaDissipation_, omega_);

    // Fix wall-cell rows to impose omega = omegaWall
    matrixConstruct_->setValues
    (
        wallCellIndices_,
        wallCellOmega_,
        wallCellFraction_
    );

    EigenVectorMap omegaSolution(omega_.data(), eIdx(numCells));

    dissipationSolver_.solve(omegaSolution, matrixA, vectorB);

    if (debug_)
    {
        const SolvePerformance& omegaPerformance =
            dissipationSolver_.lastPerformance();

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
    const Count numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar omegaLowerBound =
            k_[cellIdx] / (maxViscosityRatio_ * nu_ + vSmallValue);

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
    const Count numCells = mesh_.numCells();
    const ScalarField GammaK = computeGammaK(f1);
    cellToFaceDiffusion(GammaK, gammaKFace_);

    const ScalarField kSource{S(0.0)};

    TransportEquation equationK
    {
        .field      = Field::k,
        .phi        = k_,
        .convection = ConvectionTerm{flowRateFace, kConvectionScheme_},
        .GammaFace  = gammaKFace_,
        .source     = kSource,
        .gradPhi    = gradK_,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_->buildMatrix(equationK);

    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    // Add k source terms
    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar cellVolume = mesh_.cells()[cellIdx].volume();

        vectorB(eIdx(cellIdx)) += Pk[cellIdx] * cellVolume;

        // Destruction term: -β*·kω
        const Scalar destruction = coeffs_.betaStar * omega_[cellIdx];
        matrixA.coeffRef(eIdx(cellIdx),eIdx(cellIdx)) +=
            destruction * cellVolume;

        // -(2/3)*divU SuSp term (continuity correction)
        const Scalar suspK = (S(2.0) / S(3.0)) * divU[cellIdx];

        matrixA.coeffRef(eIdx(cellIdx),eIdx(cellIdx)) +=
            std::max(suspK, S(0.0)) * cellVolume;

        vectorB(eIdx(cellIdx)) +=
            std::max(-suspK, S(0.0)) * k_[cellIdx] * cellVolume;
    }

    // Apply implicit under-relaxation (Patankar's method)
    matrixConstruct_->relax(alphaK_, k_);

    EigenVectorMap kSolution(k_.data(), eIdx(numCells));

    kSolver_.solve(kSolution, matrixA, vectorB);

    if (debug_)
    {
        const SolvePerformance& kPerformance = kSolver_.lastPerformance();

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
    const Count numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        k_[cellIdx] = std::max(k_[cellIdx], smallValue);
    }
}


ScalarField kOmegaSST::computeTurbulentViscosity
(
    const ScalarField& f23,
    const ScalarField& strainRateMag
) const
{
    // SST turbulent viscosity:
    // nut = a1*k / max(a1*omega, b1*F23*sqrt(S2))
    const Count numCells = mesh_.numCells();
    ScalarField nut;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nut[cellIdx] =
            (coeffs_.a1 * k_[cellIdx])
          / std::max
            (
                coeffs_.a1 * omega_[cellIdx],
                f23[cellIdx] * strainRateMag[cellIdx]
            );
    }

    return nut;
}


void kOmegaSST::logFieldDiagnostics() const
{
    if (!debug_) return;

    const Count numCells = mesh_.numCells();
    if (numCells == 0) return;

    Scalar kMin = k_[0];
    Scalar kMax = k_[0];
    Scalar kSum = S(0.0);

    Scalar omegaMin = omega_[0];
    Scalar omegaMax = omega_[0];
    Scalar omegaSum = S(0.0);

    Scalar nutMin = nut_[0];
    Scalar nutMax = nut_[0];
    Scalar nutSum = S(0.0);

    #pragma omp parallel for schedule(static) \
        reduction(+:kSum, omegaSum, nutSum) \
        reduction(min:kMin, omegaMin, nutMin) \
        reduction(max:kMax, omegaMax, nutMax)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        kMin = std::min(kMin, k_[cellIdx]);
        kMax = std::max(kMax, k_[cellIdx]);
        kSum += k_[cellIdx];

        omegaMin = std::min(omegaMin, omega_[cellIdx]);
        omegaMax = std::max(omegaMax, omega_[cellIdx]);
        omegaSum += omega_[cellIdx];

        nutMin = std::min(nutMin, nut_[cellIdx]);
        nutMax = std::max(nutMax, nut_[cellIdx]);
        nutSum += nut_[cellIdx];
    }

    const Scalar n = S(numCells);

    Logger::subsection("Turbulence field statistics");
    Logger::scalarStat("k", kMin, kMax, kSum / n);
    Logger::scalarStat("omega", omegaMin, omegaMax, omegaSum / n);
    Logger::scalarStat("nut", nutMin, nutMax, nutSum / n);
}
