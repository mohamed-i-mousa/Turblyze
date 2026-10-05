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
#include <cmath>
#include <format>
#include <iostream>

// Project headers
#include "Scalar.h"
#include "HaloExchange.h"
#include "Reduce.h"
#include "Logger.h"
#include "TimeScheme.h"
#include "TurbulenceModel.h"
#include "RuntimeSelection.h"
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
    const Vector& initialVelocity,
    Scalar initialPressure,
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
    Ux_.setAll(initialVelocity.x());
    Uy_.setAll(initialVelocity.y());
    Uz_.setAll(initialVelocity.z());
    p_.setAll(initialPressure);

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
    const Vector& initialVelocity,
    Scalar initialPressure,
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
                initialVelocity,
                initialPressure,
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
                initialVelocity,
                initialPressure,
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
        Logger::sectionHeader("Starting " + algorithmName() + " Loop");
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
