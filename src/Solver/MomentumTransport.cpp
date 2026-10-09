/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MomentumTransport.cpp
 * @brief Shared methods for pressure-velocity coupling algorithms
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "MomentumTransport.h"

// Standard library headers
#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>

// Project headers
#include "Scalar.h"
#include "HaloExchange.h"
#include "Reduce.h"
#include "Logger.h"
#include "LinearInterpolation.h"
#include "TimeScheme.h"
#include "TurbulenceModel.h"
#include "RuntimeSelection.h"
#include "Initializer.h"
#include "SIMPLE.h"
#include "PISO.h"

// ************************* Special Member Functions *************************

MomentumTransport::MomentumTransport
(
    const Mesh& mesh,
    BoundaryConditions& bc,
    const TimeScheme& timeScheme,
    const GradientScheme& gradScheme,
    TurbulenceModel& turbulence,
    const Initializer& initializer,
    Scalar deltaT,
    Scalar rho,
    Scalar mu,
    Count maxIterations,
    Scalar convergenceTolerance,
    Count nOuterCorrectors,
    bool debug
)
:
    mesh_{mesh},
    bcManager_{bc},
    timeScheme_{timeScheme},
    gradientScheme_{gradScheme},
    turbulence_{turbulence},
    nu_{mu / rho},
    deltaT_{deltaT},
    maxIterations_{maxIterations},
    nOuterCorrectors_{nOuterCorrectors},
    tolerance_{convergenceTolerance},
    debug_{debug}
{
    initializer.initialize(mesh_, bcManager_, Ux_, Uy_, Uz_, p_);
    updateSymmetryBoundaries();

    // Initialize RhieChowFlowRate_ with linear interpolation
    const Count numFaces = mesh_.numFaces();

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];
        Vector Uf;

        if (face.isBoundary())
        {
            const Index owner = face.ownerCell();
            Uf = Vector
            (
                bcManager_.faceValue(face, Ux_[owner], Field::Ux),
                bcManager_.faceValue(face, Uy_[owner], Field::Uy),
                bcManager_.faceValue(face, Uz_[owner], Field::Uz)
            );
            RhieChowFlowRate_[faceIdx] =
                (face.patch()->type() != PatchType::processor
              && bcManager_.boundaryType(face, Field::Ux).isSymmetry())
              ? S(0.0)
              : dot(Uf, face.normal() * face.projectedArea());
        }
        else
        {
            Uf = Vector
            (
                interpolateToFace(mesh_, face, Ux_),
                interpolateToFace(mesh_, face, Uy_),
                interpolateToFace(mesh_, face, Uz_)
            );

            const Vector Sf = face.normal() * face.projectedArea();
            RhieChowFlowRate_[faceIdx] = dot(Uf, Sf);
        }

        UxAvgf_[faceIdx] = Uf.x();
        UyAvgf_[faceIdx] = Uf.y();
        UzAvgf_[faceIdx] = Uf.z();
    }

    // Collective: constructed on every rank together
    totalDomainCells_ = globalSum(mesh_.numDomainCells());
}

// **************************** Runtime Selection *****************************

std::unique_ptr<MomentumTransport> MomentumTransport::create
(
    const Name& algorithm,
    const Mesh& mesh,
    BoundaryConditions& bc,
    const TimeScheme& timeScheme,
    const GradientScheme& gradScheme,
    const ConvectionScheme& momentumConvectionScheme,
    LinearSolver& momentumSolver,
    LinearSolver& pressureSolver,
    TurbulenceModel& turbulence,
    const Initializer& initializer,
    Scalar deltaT,
    Scalar rho,
    Scalar mu,
    Scalar alphaU,
    Scalar alphaP,
    Count maxIterations,
    Scalar convergenceTolerance,
    Count nNonOrthogonalCorrectors,
    Count nOuterCorrectors,
    Count nPrimeCorrectors,
    bool debug
)
{
    if (algorithm == "SIMPLE")
    {
        return
            std::make_unique<SIMPLE>
            (
                mesh,
                bc,
                timeScheme,
                gradScheme,
                momentumConvectionScheme,
                momentumSolver,
                pressureSolver,
                turbulence,
                initializer,
                deltaT,
                rho,
                mu,
                alphaU,
                alphaP,
                maxIterations,
                convergenceTolerance,
                nNonOrthogonalCorrectors,
                nOuterCorrectors,
                debug
            );
    }

    if (algorithm == "PISO")
    {
        return
            std::make_unique<PISO>
            (
                mesh,
                bc,
                timeScheme,
                gradScheme,
                momentumConvectionScheme,
                momentumSolver,
                pressureSolver,
                turbulence,
                initializer,
                deltaT,
                rho,
                mu,
                alphaU,
                alphaP,
                maxIterations,
                convergenceTolerance,
                nNonOrthogonalCorrectors,
                nOuterCorrectors,
                nPrimeCorrectors,
                debug
            );
    }

    RuntimeSelection::unknownSelection
    (
        "solution algorithm",
        algorithm,
        availableAlgorithms()
    );
}


NameList MomentumTransport::availableAlgorithms()
{
    return {"SIMPLE", "PISO"};
}

// *********************************** Solve **********************************

void MomentumTransport::solve
(
    Count step,
    Count totalSteps,
    Scalar time,
    TransientFields* prevStep
)
{
    if (prevStep)
    {
        // Previous-step fields phi^n for the transient term
        prevStep->UxPrevStep = Ux_;
        prevStep->UyPrevStep = Uy_;
        prevStep->UzPrevStep = Uz_;

        // Converged t^n face flux for the Rhie-Chow
        prevStep->fluxPrevStep = faceMassFlux();

        // Previous-time-step turbulence fields
        turbulence_.beginTimeStep();
    }
    else
    {
        Logger::sectionHeader(std::format("Starting {} Loop", algorithmName()));
    }

    reportPerIteration_ = prevStep ? debug_ : true;

    const Count maxIters = prevStep ? nOuterCorrectors_ : maxIterations_;

    // Reset first-iteration residual references
    massImbalance0_ = S(0.0);
    velocityResidual0_ = S(0.0);
    pressureResidual0_ = S(0.0);
    turbulenceResidual0_.clear();

    Count iteration = 0;
    bool converged = false;

    while (iteration < maxIters && !converged)
    {
        if (debug_)
        {
            Logger::iterationHeader(iteration + 1);
            Logger::residualTableHeader();
        }
        else if (reportPerIteration_)
        {
            std::cout << std::format(" Iteration {}\n", iteration + 1);
        }

        // Previous-iteration velocity for the velocity residual
        UxPrevIter_ = Ux_;
        UyPrevIter_ = Uy_;
        UzPrevIter_ = Uz_;

        converged = outerIteration(prevStep);

        if (debug_)
        {
            Logger::iterationFooter();
        }

        iteration++;
    }

    if (prevStep)
    {
        // Roll the Crank-Nicolson stored time derivatives forward one step
        updatePrevStepDerivatives(*prevStep);
        turbulence_.updatePrevStepDerivatives();

        const CourantNumber courant = computeCourant();

        std::cout << std::format
        (
            " Time = {:.3e} s   step {}/{}   Courant max = {:.3e} mean = {:.3e}\n",
            time, step, totalSteps, courant.max, courant.mean
        );

        Logger::residualSummary
        (
            lastScaledMass_,
            lastScaledVelocity_,
            lastScaledPressure_,
            lastScaledTurbulence_
        );
    }
    else if (converged)
    {
        std::cout << std::format
        (
            "{} algorithm converged in {} iterations.\n",
            algorithmName(), iteration
        );
    }
    else
    {
        std::cout << std::format
        (
            "WARNING: {} algorithm did not converge after {} iterations.\n",
            algorithmName(), maxIterations_
        );
    }
}


bool MomentumTransport::isTransient() const noexcept
{
    return timeScheme_.isTransient();
}

// ****************************** Shared Helpers ******************************

void MomentumTransport::updatePrevStepDerivatives(TransientFields& prevStep)
{
    // Called only on the transient path
    const Count numCells = mesh_.numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar volume = mesh_.cells()[cellIdx].volume();

        prevStep.UxDdtPrevStep[cellIdx] = 
            timeScheme_.updateDdtPrevStep
            (
                volume,
                deltaT_,
                Ux_[cellIdx],
                prevStep.UxPrevStep[cellIdx],
                prevStep.UxDdtPrevStep[cellIdx]
            );

        prevStep.UyDdtPrevStep[cellIdx] =
            timeScheme_.updateDdtPrevStep
            (
                volume,
                deltaT_,
                Uy_[cellIdx],
                prevStep.UyPrevStep[cellIdx],
                prevStep.UyDdtPrevStep[cellIdx]
            );

        prevStep.UzDdtPrevStep[cellIdx] =
            timeScheme_.updateDdtPrevStep
            (
                volume,
                deltaT_,
                Uz_[cellIdx],
                prevStep.UzPrevStep[cellIdx],
                prevStep.UzDdtPrevStep[cellIdx]
            );
    }
}


void MomentumTransport::updateSymmetryBoundaries()
{
    bcManager_.refresh(mesh_, Ux_, Uy_, Uz_);
}


void MomentumTransport::updateVelocityGradients()
{
    updateSymmetryBoundaries();

    const Count numDomainCells = mesh_.numDomainCells();

    for (Index cellIdx = 0; cellIdx < numDomainCells; ++cellIdx)
    {
        gradUx_[cellIdx] =
            gradientScheme_.cellGradient(Field::Ux, Ux_, cellIdx);
        gradUy_[cellIdx] =
            gradientScheme_.cellGradient(Field::Uy, Uy_, cellIdx);
        gradUz_[cellIdx] =
            gradientScheme_.cellGradient(Field::Uz, Uz_, cellIdx);
    }

    gradientScheme_.limitGradient(Field::Ux, Ux_, gradUx_);
    gradientScheme_.limitGradient(Field::Uy, Uy_, gradUy_);
    gradientScheme_.limitGradient(Field::Uz, Uz_, gradUz_);

    // Deferred correction reads both cells of every cut face
    Halo::exchange({&gradUx_, &gradUy_, &gradUz_});

    // Assembling exchanged components replaces an exchange
    const Count numCells = mesh_.numCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradU_[cellIdx] =
            tensorFromRows
            (
                gradUx_[cellIdx],
                gradUy_[cellIdx],
                gradUz_[cellIdx]
            );
    }
}


void MomentumTransport::solveTurbulence()
{
    if (!turbulence_.isTurbulent())
    {
        return;
    }

    updateVelocityGradients();

    turbulence_.solve
    (
        Ux_,
        Uy_,
        Uz_,
        faceMassFlux(),
        gradU_
    );
}


bool MomentumTransport::checkConvergence()
{
    const TurbulenceModel::ResidualPair turbulenceResiduals =
        turbulence_.residualOutputs();

    // Compute raw residuals
    const Scalar massImbalance = this->massImbalance();
    const Scalar velocityResidual = this->velocityResidual();
    const Scalar pressureResidual = this->pressureResidual();

    Scalar scaledMass = massImbalance;
    Scalar scaledVelocity = velocityResidual;
    Scalar scaledPressure = pressureResidual;
    std::vector<Logger::Residuals> scaledTurbulenceResiduals;
    bool converged = false;

    if (isTransient())
    {
        // In transient simulations, report the inherent dimensionless residuals directly.
        converged =
            (scaledMass < tolerance_)
         && (scaledVelocity < tolerance_)
         && (scaledPressure < tolerance_);

        if (turbulence_.isTurbulent())
        {
            scaledTurbulenceResiduals.reserve(turbulenceResiduals.size());
            for (const auto& residual : turbulenceResiduals)
            {
                scaledTurbulenceResiduals.push_back
                (
                    {residual.first, residual.second}
                );
                converged = converged && (residual.second < tolerance_);
            }
        }
    }
    else
    {
        // Steady-state simulations (e.g. SIMPLE): scale relative to the initial iteration
        if (massImbalance0_ < vSmallValue)
        {
            massImbalance0_ = massImbalance;
            velocityResidual0_ = velocityResidual;
            pressureResidual0_ = pressureResidual;

            turbulenceResidual0_.clear();
            for (const auto& residual : turbulenceResiduals)
            {
                turbulenceResidual0_.push_back(residual.second);
            }
        }

        scaledMass = massImbalance / (massImbalance0_ + vSmallValue);
        scaledVelocity = velocityResidual / (velocityResidual0_ + vSmallValue);
        scaledPressure = pressureResidual / (pressureResidual0_ + vSmallValue);

        converged =
            (scaledMass < tolerance_)
         && (scaledVelocity < tolerance_)
         && (scaledPressure < tolerance_);

        if (turbulence_.isTurbulent())
        {
            const Count residualCount =
                std::min
                (
                    turbulenceResiduals.size(),
                    turbulenceResidual0_.size()
                );

            scaledTurbulenceResiduals.reserve(residualCount);

            for (Index i = 0; i < residualCount; ++i)
            {
                const Scalar scaled =
                    turbulenceResiduals[i].second
                  / (turbulenceResidual0_[i] + vSmallValue);

                scaledTurbulenceResiduals.push_back
                (
                    {turbulenceResiduals[i].first, scaled}
                );

                converged = converged && (scaled < tolerance_);
            }

            if (residualCount != turbulenceResiduals.size())
            {
                converged = false;
            }
        }
    }

    // Remember the latest scaled residuals for the per-time-step summary
    lastScaledMass_ = scaledMass;
    lastScaledVelocity_ = scaledVelocity;
    lastScaledPressure_ = scaledPressure;
    lastScaledTurbulence_ = scaledTurbulenceResiduals;

    if (debug_)
    {
        Logger::subsection(isTransient() ? "Residuals" : "Scaled residuals");
        Logger::scaledResidual("mass",     scaledMass);
        Logger::scaledResidual("velocity", scaledVelocity);
        Logger::scaledResidual("pressure", scaledPressure);
        if (turbulence_.isTurbulent())
        {
            for (const Logger::Residuals& residual : scaledTurbulenceResiduals)
            {
                Logger::scaledResidual(residual.first, residual.second);
            }
        }
    }
    else if (reportPerIteration_)
    {
        // Laminar runs carry an empty turbulence span, printing nothing extra
        Logger::residualSummary
        (
            scaledMass,
            scaledVelocity,
            scaledPressure,
            scaledTurbulenceResiduals
        );
    }

    return converged;
}


Scalar MomentumTransport::massImbalance() const noexcept
{
    // Dimensionless normalized continuity residual per cell, averaged
    const FaceFluxField& flux = faceMassFlux();

    Scalar totalNormImbalance = S(0.0);

    const Count numCells = mesh_.numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const auto& cell = mesh_.cells()[cellIdx];
        const auto& faceIndices = cell.faceIndices();
        const auto& faceSigns = cell.faceSigns();

        Scalar net = S(0.0);
        Scalar sumAbs = S(0.0);

        for (Index j = 0; j < faceIndices.size(); ++j)
        {
            const Index faceIdx = faceIndices[j];
            const int sign = faceSigns[j];
            const Scalar mf = flux[faceIdx];
            net += S(sign) * mf;
            sumAbs += std::abs(mf);
        }

        const Scalar denom = sumAbs + vSmallValue;
        totalNormImbalance += std::abs(net) / denom;
    }

    return globalSum(totalNormImbalance)
         / S(std::max<Count>(1, totalDomainCells_));
}


Scalar MomentumTransport::velocityResidual() const noexcept
{
    // Normalized residual: ||U - U_prev||_2 / (||U_prev||_2 + eps)
    Scalar num = S(0.0);
    Scalar den = S(0.0);

    const Count numCells = mesh_.numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar dx = Ux_[cellIdx] - UxPrevIter_[cellIdx];
        const Scalar dy = Uy_[cellIdx] - UyPrevIter_[cellIdx];
        const Scalar dz = Uz_[cellIdx] - UzPrevIter_[cellIdx];

        num += dx * dx + dy * dy + dz * dz;
        den += UxPrevIter_[cellIdx] * UxPrevIter_[cellIdx]
             + UyPrevIter_[cellIdx] * UyPrevIter_[cellIdx]
             + UzPrevIter_[cellIdx] * UzPrevIter_[cellIdx];
    }

    num = std::sqrt(globalSum(num) + vSmallValue);
    den = std::sqrt(globalSum(den) + vSmallValue);

    return num / den;
}


MomentumTransport::CourantNumber
MomentumTransport::computeCourant() const noexcept
{
    const FaceFluxField& flux = faceMassFlux();

    const Count numCells = mesh_.numDomainCells();

    Scalar maxCourant = S(0.0);
    Scalar sumCourant = S(0.0);

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const auto& cell = mesh_.cells()[cellIdx];
        const auto& faceIndices = cell.faceIndices();

        Scalar sumFlux = S(0.0);
        for (Index j = 0; j < faceIndices.size(); ++j)
        {
            sumFlux += std::abs(flux[faceIndices[j]]);
        }

        const Scalar courant =
            S(0.5) * sumFlux * deltaT_ / cell.volume();
        maxCourant = std::max(maxCourant, courant);
        sumCourant += courant;
    }

    return
    {
        globalMax(maxCourant),
        globalSum(sumCourant) / S(std::max<Count>(1, totalDomainCells_))
    };
}


// ********************** Shared Finite-Volume Methods ***********************

void MomentumTransport::updateEffectiveViscosity()
{
    const Count numCells = mesh_.numCells();
    const Count numFaces = mesh_.numFaces();

    const ScalarField& nut = turbulence_.turbulentViscosity();

    // Build cell-based effective viscosity
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nuEff_[cellIdx] = nu_ + nut[cellIdx];
    }

    // Build face-based effective viscosity
    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            // Turbulent models may provide wall-function boundary nut.
            nuEffFace_[faceIdx] =
                nu_
              + turbulence_.boundaryTurbulentViscosity(face);
        }
        else
        {
            // Internal faces: linear interpolation
            nuEffFace_[faceIdx] = interpolateToFace(mesh_, face, nuEff_);
        }
    }
}


void MomentumTransport::buildFaceDiagonal()
{
    const Count numFaces = mesh_.numFaces();

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            // Dirichlet p' couples pressure and velocity through the face; a
            // zero-gradient or symmetry plane decouples them
            DUf_[faceIdx] =
                (face.patch()->type() != PatchType::processor
              && bcManager_.boundaryType(face.patch()->name(), Field::pCorr).fixesValue())
              ? DU_[face.ownerCell()]
              : S(0.0);
        }
        else
        {
            // Internal faces
            DUf_[faceIdx] = interpolateToFace(mesh_, face, DU_);
        }
    }
}


void MomentumTransport::updateRhieChowFlowRate
(
    Scalar alphaU,
    const TransientFields* prevStep
)
{
    const Count numFaces = mesh_.numFaces();

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            const Index owner = face.ownerCell();
            const Vector Uf
            (
                bcManager_.faceValue(face, Ux_[owner], Field::Ux),
                bcManager_.faceValue(face, Uy_[owner], Field::Uy),
                bcManager_.faceValue(face, Uz_[owner], Field::Uz)
            );
            UxAvgf_[faceIdx] = Uf.x();
            UyAvgf_[faceIdx] = Uf.y();
            UzAvgf_[faceIdx] = Uf.z();

            RhieChowFlowRate_[faceIdx] =
                (face.patch()->type() != PatchType::processor
              && bcManager_.boundaryType(face, Field::Ux).isSymmetry())
              ? S(0.0)
              : dot(Uf, face.normal() * face.projectedArea());
            continue;
        }

        const Index P = face.ownerCell();
        const Index N = face.neighborCell().value();

        // Linear-interpolated velocity at face
        const Vector UfLinear
        (
            interpolateToFace(mesh_, face, Ux_),
            interpolateToFace(mesh_, face, Uy_),
            interpolateToFace(mesh_, face, Uz_)
        );

        const Vector gradPAvgf = interpolateToFace(mesh_, face, gradP_);
        const Vector Sf = face.normal() * face.projectedArea();
        const Vector gradPf =
            gradientScheme_.faceGradient
            (
                p_,
                gradP_[P],
                gradP_[N],
                faceIdx
            );
        const Vector UfPrevIter
        (
            UxAvgPrevIterf_[faceIdx],
            UyAvgPrevIterf_[faceIdx],
            UzAvgPrevIterf_[faceIdx]
        );

        RhieChowFlowRate_[faceIdx] =
            dot(UfLinear, Sf)
          - dot((DUf_[faceIdx] * (gradPf - gradPAvgf)), Sf)
          + (S(1.0) - alphaU)
          * (RhieChowFlowRatePrevIter_[faceIdx] - dot(UfPrevIter, Sf));

        // prevStep is non-null exactly on the transient path
        if (prevStep != nullptr)
        {
            const Vector UfPrevStepLinear
            (
                interpolateToFace(mesh_, face, prevStep->UxPrevStep),
                interpolateToFace(mesh_, face, prevStep->UyPrevStep),
                interpolateToFace(mesh_, face, prevStep->UzPrevStep)
            );
            const Scalar phiCorr =
                prevStep->fluxPrevStep[faceIdx] - dot(UfPrevStepLinear, Sf);
            const Scalar coeff = S(1.0) - std::min
            (
                std::abs(phiCorr)
              / (std::abs(prevStep->fluxPrevStep[faceIdx]) + vSmallValue),
                S(1.0)
            );
            const Scalar DTf = DUf_[faceIdx] * coeff / deltaT_;
            RhieChowFlowRate_[faceIdx] += DTf * phiCorr;
        }
    }
}


void MomentumTransport::addTransposeGradientSource
(
    ScalarField& UxSource,
    ScalarField& UySource,
    ScalarField& UzSource
) const
{
    const Count numCells = mesh_.numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar sumX = S(0.0);
        Scalar sumY = S(0.0);
        Scalar sumZ = S(0.0);

        const auto& cell = mesh_.cells()[cellIdx];
        const auto& faceIndices = cell.faceIndices();
        const auto& faceSigns = cell.faceSigns();

        for (Index j = 0; j < faceIndices.size(); ++j)
        {
            const Index faceIdx = faceIndices[j];
            const Face& face = mesh_.faces()[faceIdx];

            // Symmetry planes integrate full viscous normal stress implicitly;
            // skip to avoid double-counting the transpose stress here
            if
            (
                face.isBoundary()
             && face.patch()->type() != PatchType::processor
             && bcManager_.boundaryType(face, Field::Ux).isSymmetry()
            )
            {
                continue;
            }

            const Scalar sign = S(faceSigns[j]);
            const Vector Sf = face.normal() * face.projectedArea() * sign;
            const Scalar nuEfff = nuEffFace_[faceIdx];

            Tensor gradUf;
            if (face.isBoundary())
            {
                gradUf = gradU_[cellIdx];
            }
            else
            {
                gradUf = interpolateToFace(mesh_, face, gradU_);
            }

            sumX += nuEfff * dot(gradUf.col(0), Sf);
            sumY += nuEfff * dot(gradUf.col(1), Sf);
            sumZ += nuEfff * dot(gradUf.col(2), Sf);
        }

        UxSource[cellIdx] += sumX;
        UySource[cellIdx] += sumY;
        UzSource[cellIdx] += sumZ;
    }
}


