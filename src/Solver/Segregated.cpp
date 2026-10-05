/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Segregated.cpp
 * @brief Pressure-correction and Rhie-Chow methods for segregated algorithms
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Segregated.h"

// Standard library headers
#include <cmath>
#include <iostream>
#include <algorithm>

// Project headers
#include "Scalar.h"
#include "HaloExchange.h"
#include "PETScRuntime.h"
#include "Reduce.h"
#include "Logger.h"
#include "LinearInterpolation.h"
#include "TimeScheme.h"
#include "TurbulenceModel.h"

// ************************* Special Member Functions *************************

Segregated::Segregated
(
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
    bool debug
)
:
    MomentumTransport
    {
        mesh,
        bc,
        timeScheme,
        gradScheme,
        turbulence,
        initialVelocity,
        initialPressure,
        deltaT,
        rho,
        mu,
        maxIterations,
        convergenceTolerance,
        nOuterCorrectors,
        debug
    },
    momentumConvectionScheme_{momentumConvectionScheme},
    momentumSolver_{momentumSolver},
    pressureSolver_{pressureSolver},
    matrixConstruct_{mesh, bc},
    alphaU_{alphaU},
    alphaP_{alphaP},
    nNonOrthogonalCorrectors_{nNonOrthogonalCorrectors}
{
    // No Dirichlet anchor on ANY rank leaves p' with a constant null space
    Count fixedPressurePatches = 0;

    for (const BoundaryPatch& patch : mesh.patches())
    {
        if (patch.type() == PatchType::processor)
        {
            continue;
        }

        if (bc.boundaryType(patch.name(), Field::p).fixesValue())
        {
            ++fixedPressurePatches;
        }
    }

    pCorrNeedsNullSpace_ = globalSum(fixedPressurePatches) == 0;

    // Update boundary coefficients
    updateSymmetryBoundaries();

    UxAvgf_.setAll(initialVelocity.x());
    UyAvgf_.setAll(initialVelocity.y());
    UzAvgf_.setAll(initialVelocity.z());

    // Initialize RhieChowFlowRate_ with linear interpolation
    const Count numFaces = mesh.numFaces();

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh.faces()[faceIdx];
        Vector Uf;

        if (face.isBoundary())
        {
            const Index owner = face.ownerCell();
            Uf = Vector
            (
                bcManager().faceValue(face, Ux()[owner], Field::Ux),
                bcManager().faceValue(face, Uy()[owner], Field::Uy),
                bcManager().faceValue(face, Uz()[owner], Field::Uz)
            );
        }
        else
        {
            Uf = Vector
            (
                interpolateToFace(mesh, face, Ux()),
                interpolateToFace(mesh, face, Uy()),
                interpolateToFace(mesh, face, Uz())
            );
        }

        const bool isZeroFlux =
            face.isBoundary()
         && face.patch()->type() != PatchType::processor
         && bcManager().boundaryType(face.patch()->name(), Field::Ux).isSymmetry();

        const Vector Sf = face.normal() * face.projectedArea();
        RhieChowFlowRate_[faceIdx] = isZeroFlux ? S(0.0) : dot(Uf, Sf);
    }
}


Scalar Segregated::pressureResidual() const noexcept
{
    // Normalize p' RMS by RMS(p)
    Scalar sumP2 = S(0.0);

    const Count numCells = mesh().numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        sumP2 += pressure()[cellIdx] * pressure()[cellIdx];
    }

    const Scalar pRms =
        std::sqrt(globalSum(sumP2) / S(totalDomainCells()));

    return lastPressureCorrectionRMS_ / (pRms + vSmallValue);
}


void Segregated::updateEffectiveViscosity()
{
    const Count numCells = mesh().numCells();
    const Count numFaces = mesh().numFaces();

    const ScalarField& nut = turbulence().turbulentViscosity();

    // Build cell-based effective viscosity
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        nuEff_[cellIdx] = nu() + nut[cellIdx];
    }

    // Build face-based effective viscosity
    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh().faces()[faceIdx];

        if (face.isBoundary())
        {
            // Turbulent models may provide wall-function boundary nut.
            nuEffFace_[faceIdx] =
                nu()
              + turbulence().boundaryTurbulentViscosity(face);
        }
        else
        {
            // Internal faces: linear interpolation
            nuEffFace_[faceIdx] = interpolateToFace(mesh(), face, nuEff_);
        }
    }
}


void Segregated::assembleMomentum()
{
    const Count numCells = mesh().numDomainCells();

    // Reset diagonals accumulator
    DU_.setAll(S(0.0));

    // Pressure gradient source term (depends on the current gradP_)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar volume = mesh().cells()[cellIdx].volume();
        UxSource_[cellIdx] = -gradP_[cellIdx].x() * volume;
        UySource_[cellIdx] = -gradP_[cellIdx].y() * volume;
        UzSource_[cellIdx] = -gradP_[cellIdx].z() * volume;
    }

    updateVelocityGradients();

    // Add transpose gradient source term (only varies with turbulent nut)
    if (turbulence().isTurbulent())
    {
        addTransposeGradientSource();
    }
}


void Segregated::solveMomentum(const TransientFields* prevStep)
{
    // Cache face velocities and mass flux for the next assembly
    UxAvgPrevIterf_ = UxAvgf_;
    UyAvgPrevIterf_ = UyAvgf_;
    UzAvgPrevIterf_ = UzAvgf_;
    RhieChowFlowRatePrevIter_ = RhieChowFlowRate_;

    // Seed the pressure gradient for the momentum source and Rhie-Chow
    gradientScheme().fieldGradient(Field::p, pressure(), gradP_);
    Halo::exchange({&gradP_});

    updateEffectiveViscosity();
    assembleMomentum();

    const ConvectionTerm convection
    {
        RhieChowFlowRatePrevIter_,
        momentumConvectionScheme_
    };


    TransportEquation equations[]
    {
        {
            .field          = Field::Ux,
            .phi            = Ux(),
            .transient      = ddtTerm(prevStep,
                &TransientFields::UxPrevStep, &TransientFields::UxDdtPrevStep),
            .convection     = convection,
            .GammaFace      = nuEffFace_,
            .source         = UxSource_,
            .gradPhi        = gradUx(),
            .gradScheme     = gradientScheme()
        },
        {
            .field          = Field::Uy,
            .phi            = Uy(),
            .transient      = ddtTerm(prevStep,
                &TransientFields::UyPrevStep, &TransientFields::UyDdtPrevStep),
            .convection     = convection,
            .GammaFace      = nuEffFace_,
            .source         = UySource_,
            .gradPhi        = gradUy(),
            .gradScheme     = gradientScheme()
        },
        {
            .field          = Field::Uz,
            .phi            = Uz(),
            .transient      = ddtTerm(prevStep,
                &TransientFields::UzPrevStep, &TransientFields::UzDdtPrevStep),
            .convection     = convection,
            .GammaFace      = nuEffFace_,
            .source         = UzSource_,
            .gradPhi        = gradUz(),
            .gradScheme     = gradientScheme()
        }
    };

    const ScalarField* prevIters[]
    {
        &UxPrevIter(),
        &UyPrevIter(),
        &UzPrevIter()
    };

    // Build and implicitly solve each under-relaxed component
    for
    (
        Index momentumComponent = 0;
        momentumComponent < 3;
        ++momentumComponent
    )
    {
        matrixConstruct_.buildMatrix(equations[momentumComponent]);

        matrixConstruct_.relax(alphaU_, *prevIters[momentumComponent]);

        diagonalDU(momentumComponent);

        matrixConstruct_.assemble();


        momentumSolver_.solve
        (
            {equations[momentumComponent].phi.data(), mesh().numDomainCells()},
            matrixConstruct_.matrixA(),
            matrixConstruct_.rhsVec()
        );

        if (debug())
        {
            const SolvePerformance& momentumPerformance =
                momentumSolver_.lastPerformance();

            Logger::residualRow
            (
                fieldToString(equations[momentumComponent].field),
                momentumPerformance.solverName,
                momentumPerformance.iterations,
                momentumPerformance.finalResidual
            );
        }
    }

    // KSP writes owned entries only: refresh U ghosts before any face read
    Halo::exchange({&Ux(), &Uy(), &Uz()});

    buildFaceDiagonal();
}


void Segregated::diagonalDU(Index component)
{
    const std::span<const Scalar> diagonal = matrixConstruct_.diagonal();
    const Count numCells = mesh().numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        DU_[cellIdx] += diagonal[cellIdx];
    }

    if (component == 2)
    {
        for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            DU_[cellIdx] =
                S(3.0) * mesh().cells()[cellIdx].volume()
              / (DU_[cellIdx] + vSmallValue);
        }

        // Rhie-Chow and the p' diffusion interpolate DU at cut faces
        Halo::exchange({&DU_});
    }
}


void Segregated::buildFaceDiagonal()
{
    const Count numFaces = mesh().numFaces();

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh().faces()[faceIdx];

        if (face.isBoundary())
        {
            // Dirichlet p' couples pressure and velocity through the face; a
            // zero-gradient or symmetry plane decouples them
            DUf_[faceIdx] =
                (face.patch()->type() != PatchType::processor
              && bcManager().boundaryType(face.patch()->name(), Field::pCorr).fixesValue())
              ? DU_[face.ownerCell()]
              : S(0.0);
        }
        else
        {
            // Internal faces
            DUf_[faceIdx] = interpolateToFace(mesh(), face, DU_);
        }
    }
}


void Segregated::updateRhieChowFlowRate(const TransientFields* prevStep)
{
    const Count numFaces = mesh().numFaces();

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh().faces()[faceIdx];

        if (face.isBoundary())
        {
            const Index owner = face.ownerCell();

            UxAvgf_[faceIdx] = bcManager().faceValue(face, Ux()[owner], Field::Ux);
            UyAvgf_[faceIdx] = bcManager().faceValue(face, Uy()[owner], Field::Uy);
            UzAvgf_[faceIdx] = bcManager().faceValue(face, Uz()[owner], Field::Uz);

            // A flux-constrained face (symmetry) carries zero mass flux
            const bool isZeroFlux =
                face.patch()->type() != PatchType::processor
             && bcManager().boundaryType(face.patch()->name(), Field::Ux).isSymmetry();

            const Vector Uf
            (
                UxAvgf_[faceIdx],
                UyAvgf_[faceIdx],
                UzAvgf_[faceIdx]
            );

            RhieChowFlowRate_[faceIdx] =
                isZeroFlux
              ? S(0.0)
              : dot(Uf, face.normal() * face.projectedArea());
            continue;
        }

        const Index P = face.ownerCell();
        const Index N = face.neighborCell().value();

        // Linear-interpolated velocity at face
        const Vector UfLinear
        (
            interpolateToFace(mesh(), face, Ux()),
            interpolateToFace(mesh(), face, Uy()),
            interpolateToFace(mesh(), face, Uz())
        );

        const Vector gradPAvgf = interpolateToFace(mesh(), face, gradP_);
        const Vector Sf = face.normal() * face.projectedArea();
        const Vector gradPf =
            gradientScheme().faceGradient
            (
                pressure(),
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
          + (S(1.0) - alphaU_)
          * (RhieChowFlowRatePrevIter_[faceIdx] - dot(UfPrevIter, Sf));

        // prevStep is non-null exactly on the transient path
        if (prevStep != nullptr)
        {
            const Vector UfPrevStepLinear
            (
                interpolateToFace(mesh(), face, prevStep->UxPrevStep),
                interpolateToFace(mesh(), face, prevStep->UyPrevStep),
                interpolateToFace(mesh(), face, prevStep->UzPrevStep)
            );
            const Scalar phiCorr =
                prevStep->fluxPrevStep[faceIdx] - dot(UfPrevStepLinear, Sf);
            const Scalar coeff = S(1.0) - std::min
            (
                std::abs(phiCorr)
              / (std::abs(prevStep->fluxPrevStep[faceIdx]) + vSmallValue),
                S(1.0)
            );
            const Scalar DTf = DUf_[faceIdx] * coeff / deltaT();
            RhieChowFlowRate_[faceIdx] += DTf * phiCorr;
        }
    }
}


void Segregated::solvePressureCorrection()
{
    const Count numCells = mesh().numDomainCells();

    // Compute mass imbalance source term
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar net = S(0.0);
        const auto& faceIndices = mesh().cells()[cellIdx].faceIndices();
        const auto& signs = mesh().cells()[cellIdx].faceSigns();

        for (Index j = 0; j < faceIndices.size(); ++j)
        {
            net += signs[j] * RhieChowFlowRate_[faceIndices[j]];
        }

        massImbalanceSrc_[cellIdx] = -net;
    }

    // p' restarts from zero every outer iteration
    pCorr_.setAll(S(0.0));
    gradPCorr_.setAll(Vector{S(0.0), S(0.0), S(0.0)});

    TransportEquation equationPCorr
    {
        .field      = Field::pCorr,
        .phi        = pCorr_,
        .convection = std::nullopt,
        .GammaFace  = DUf_,
        .source     = massImbalanceSrc_,
        .gradPhi    = gradPCorr_,
        .gradScheme = gradientScheme()
    };

    // Attached around the corrector loop only (the Matrix is shared)
    if (pCorrNeedsNullSpace_)
    {
        MatNullSpace constantNullSpace = nullptr;
        CheckPETSc
        (
            MatNullSpaceCreate
            (
                PETScRuntime::comm(),
                PETSC_TRUE,
                0,
                nullptr,
                &constantNullSpace
            )
        );
        CheckPETSc
        (
            MatSetNullSpace(matrixConstruct_.matrixA(), constantNullSpace)
        );
        CheckPETSc(MatNullSpaceDestroy(&constantNullSpace));
    }

    for
    (
        Count corrector = 0;
        corrector <= nNonOrthogonalCorrectors_;
        ++corrector
    )
    {
        matrixConstruct_.buildMatrix(equationPCorr);

        matrixConstruct_.assemble();

        pressureSolver_.solve
        (
            {pCorr_.data(), numCells},
            matrixConstruct_.matrixA(),
            matrixConstruct_.rhsVec()
        );

        if (debug())
        {
            const SolvePerformance& pressurePerformance =
                pressureSolver_.lastPerformance();

            Logger::residualRow
            (
                "p'",
                pressurePerformance.solverName,
                pressurePerformance.iterations,
                pressurePerformance.finalResidual
            );
        }

        // The corrector reads p' and grad p' at both cells of every cut
        Halo::exchange({&pCorr_});

        // grad(p') feeds the next corrector's non-orthogonal term
        for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            gradPCorr_[cellIdx] =
                gradientScheme().cellGradient(Field::pCorr, pCorr_, cellIdx);
        }

        Halo::exchange({&gradPCorr_});
    }

    if (pCorrNeedsNullSpace_)
    {
        CheckPETSc(MatSetNullSpace(matrixConstruct_.matrixA(), nullptr));
    }
}


void Segregated::correctVelocity()
{
    const Count numCells = mesh().numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Ux()[cellIdx] -= DU_[cellIdx] * gradPCorr_[cellIdx].x();
        Uy()[cellIdx] -= DU_[cellIdx] * gradPCorr_[cellIdx].y();
        Uz()[cellIdx] -= DU_[cellIdx] * gradPCorr_[cellIdx].z();
    }

    // The face averages below read the corrected U at both cells
    Halo::exchange({&Ux(), &Uy(), &Uz()});

    updateSymmetryBoundaries();

    // Update face velocities
    const Count numFaces = mesh().numFaces();

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh().faces()[faceIdx];

        if (face.isBoundary())
        {
            const Index owner = face.ownerCell();
            UxAvgf_[faceIdx] = bcManager().faceValue(face, Ux()[owner], Field::Ux);
            UyAvgf_[faceIdx] = bcManager().faceValue(face, Uy()[owner], Field::Uy);
            UzAvgf_[faceIdx] = bcManager().faceValue(face, Uz()[owner], Field::Uz);
        }
        else
        {
            UxAvgf_[faceIdx] = interpolateToFace(mesh(), face, Ux());
            UyAvgf_[faceIdx] = interpolateToFace(mesh(), face, Uy());
            UzAvgf_[faceIdx] = interpolateToFace(mesh(), face, Uz());
        }
    }
}


void Segregated::correctFlowRate()
{
    // Update mass flux on faces
    const Count numFaces = mesh().numFaces();

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh().faces()[faceIdx];

        if (face.isBoundary())
        {
            if
            (
                face.patch()->type() == PatchType::processor
             || !bcManager().boundaryType
                (
                    face.patch()->name(), Field::p
                ).correctsBoundaryFlux()
            )
            {
                continue;
            }

            const Scalar gradn =
                dot(gradPCorr_[face.ownerCell()], face.normal());

            const Scalar flowRateCorrection =
                DU_[face.ownerCell()] * gradn * face.projectedArea();

            RhieChowFlowRate_[faceIdx] -= flowRateCorrection;
            continue;
        }

        const Index ownerIdx = face.ownerCell();
        const Index neighborIdx = face.neighborCell().value();

        const Vector gradPCorrf =
            gradientScheme().faceGradient
            (
                pCorr_,
                gradPCorr_[ownerIdx],
                gradPCorr_[neighborIdx],
                faceIdx
            );

        const Vector Sf = face.normal() * face.projectedArea();
        const Scalar flowRateCorrection =
            DUf_[faceIdx] * dot(gradPCorrf, Sf);

        RhieChowFlowRate_[faceIdx] -= flowRateCorrection;
    }
}


void Segregated::correctPressure()
{
    Scalar sumSq = S(0.0);

    const Count numCells = mesh().numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        sumSq += pCorr_[cellIdx] * pCorr_[cellIdx];
    }

    lastPressureCorrectionRMS_ =
        std::sqrt(globalSum(sumSq) / S(totalDomainCells()));

    // Apply pressure correction
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        pressure()[cellIdx] += alphaP_ * pCorr_[cellIdx];
    }

    // Next iteration's gradP stencil and Rhie-Chow read p across cuts
    Halo::exchange({&pressure()});
}


void Segregated::addTransposeGradientSource()
{
    const Count numCells = mesh().numDomainCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar sumX = S(0.0);
        Scalar sumY = S(0.0);
        Scalar sumZ = S(0.0);

        const auto& cell = mesh().cells()[cellIdx];
        const auto& faceIndices = cell.faceIndices();
        const auto& faceSigns = cell.faceSigns();

        for (Index j = 0; j < faceIndices.size(); ++j)
        {
            const Index faceIdx = faceIndices[j];
            const Scalar sign = S(faceSigns[j]);
            const Face& face = mesh().faces()[faceIdx];

            const Vector Sf = face.normal() * face.projectedArea() * sign;
            const Scalar nuEfff = nuEffFace_[faceIdx];

            Tensor gradUf;
            if (face.isBoundary())
            {
                gradUf = gradU()[cellIdx];
            }
            else
            {
                gradUf = interpolateToFace(mesh(), face, gradU());
            }

            sumX += nuEfff * dot(gradUf.col(0), Sf);
            sumY += nuEfff * dot(gradUf.col(1), Sf);
            sumZ += nuEfff * dot(gradUf.col(2), Sf);
        }

        UxSource_[cellIdx] += sumX;
        UySource_[cellIdx] += sumY;
        UzSource_[cellIdx] += sumZ;
    }
}


std::optional<TransientTerm> Segregated::ddtTerm
(
    const TransientFields* prevStep,
    ScalarField TransientFields::* phiPrevStep,
    ScalarField TransientFields::* ddtPrevStep
) const
{
    if (prevStep == nullptr)
    {
        return std::nullopt;
    }

    return 
        TransientTerm
        {
            timeScheme(),
            deltaT(),
            prevStep->*phiPrevStep,
            &(prevStep->*ddtPrevStep)
        };
}
