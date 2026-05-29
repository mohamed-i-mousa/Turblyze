/******************************************************************************
 * @file SIMPLE.cpp
 * @brief Implementation of the SIMPLE algorithm for pressure-velocity coupling
 *****************************************************************************/

// ********************************** Headers *********************************

/// Implementation header
#include "SIMPLE.h"

/// Standard library headers
#include <cmath>
#include <iostream>
#include <algorithm>

/// External library headers
#include <omp.h>

/// Project headers
#include "Scalar.h"
#include "Logger.h"
#include "LinearInterpolation.h"
#include "Constraint.h"
#include "kOmegaSST.h"

// ************************* Special Member Functions *************************

SIMPLE::SIMPLE
(
    const Mesh& mesh,
    const BoundaryConditions& bc,
    const GradientScheme& gradScheme,
    const ConvectionSchemes& momentumConvectionScheme,
    LinearSolver& momentumSolver,
    LinearSolver& pressureSolver,
    kOmegaSST* turbulence,
    const Scalar rho,
    const Scalar mu,
    const Vector& initialVelocity,
    const Scalar initialPressure,
    const Scalar alphaU,
    const Scalar alphaP,
    const int maxIterations,
    const Scalar convergenceTolerance,
    const bool velocityConstraintEnabled,
    const bool pressureConstraintEnabled,
    const Scalar maxVelocityMagnitude,
    const Scalar minPressure,
    const Scalar maxPressure,
    const bool debug
)
:
    mesh_{mesh},
    bcManager_{bc},
    gradientScheme_{gradScheme},
    momentumConvectionScheme_{momentumConvectionScheme},
    momentumSolver_{momentumSolver},
    pressureSolver_{pressureSolver},
    turbulence_{turbulence},
    matrixConstruct_{mesh_, bcManager_},
    nu_{mu / rho},
    alphaU_{alphaU},
    alphaP_{alphaP},
    maxIterations_{maxIterations},
    tolerance_{convergenceTolerance},
    debug_{debug},
    constraintSystem_
    {
        Ux_,
        Uy_,
        Uz_,
        p_,
        velocityConstraintEnabled,
        pressureConstraintEnabled,
        maxVelocityMagnitude,
        minPressure,
        maxPressure
    }
{
    Ux_.setAll(initialVelocity.x());
    Uy_.setAll(initialVelocity.y());
    Uz_.setAll(initialVelocity.z());
    UxAvgf_.setAll(initialVelocity.x());
    UyAvgf_.setAll(initialVelocity.y());
    UzAvgf_.setAll(initialVelocity.z());
    p_.setAll(initialPressure);

    // Initialize RhieChowFlowRate_ with linear interpolation
    const size_t numFaces = mesh_.numFaces();

    #pragma omp parallel for schedule(static)
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];
        Vector Uf;

        if (face.isBoundary())
        {
            Uf = Vector
            (
                bcManager_.boundaryFaceValue(Field::Ux, Ux_, face),
                bcManager_.boundaryFaceValue(Field::Uy, Uy_, face),
                bcManager_.boundaryFaceValue(Field::Uz, Uz_, face)
            );
        }
        else
        {
            Uf = Vector
            (
                interpolateToFace(face, Ux_),
                interpolateToFace(face, Uy_),
                interpolateToFace(face, Uz_)
            );
        }

        const Vector Sf = face.normal() * face.projectedArea();
        RhieChowFlowRate_[faceIdx] = dot(Uf, Sf);
    }
}

// ******************************* SIMPLE Solve *******************************

void SIMPLE::solve()
{
    Logger::sectionHeader("Starting SIMPLE Loop");

    // Reset first-iteration residual references for convergence tracking
    massImbalance0_    = S(0.0);
    velocityResidual0_ = S(0.0);
    pressureResidual0_ = S(0.0);
    kResidual0_        = S(0.0);
    omegaResidual0_    = S(0.0);

    int iteration = 0;
    bool converged = false;

    while (!converged && iteration < maxIterations_)
    {
        if (debug_)
        {
            Logger::iterationHeader(iteration + 1);
            Logger::residualTableHeader();
        }
        else
        {
            std::cout << " Iteration " << iteration + 1 << '\n';
        }

        UxPrev_ = Ux_;
        UyPrev_ = Uy_;
        UzPrev_ = Uz_;
        UxAvgPrevf_ = UxAvgf_;
        UyAvgPrevf_ = UyAvgf_;
        UzAvgPrevf_ = UzAvgf_;
        RhieChowFlowRatePrev_ = RhieChowFlowRate_;

        gradientScheme_.fieldGradient(Field::p, p_, gradP_);

        solveMomentumEquations();

        updateRhieChowFlowRate();

        solvePressureCorrection();

        correctVelocity();

        correctFlowRate();

        correctPressure();

        solveTurbulence();

        converged = checkConvergence();

        if (debug_)
        {
            Logger::iterationFooter();
        }

        iteration++;
    }

    if (!converged)
    {
        std::cout
            << "WARNING: SIMPLE algorithm did not converge after "
            << maxIterations_ << " iterations." << '\n';
    }
    else
    {
        std::cout
            << "SIMPLE algorithm converged in " << iteration
            << " iterations." << '\n';
    }
}

// ****************************** Private Methods *****************************

void SIMPLE::solveMomentumEquations()
{
    const size_t numCells = mesh_.numCells();
    const size_t numFaces = mesh_.numFaces();

    // Reset diagonals accumulator
    DU_.setAll(S(0.0));
    DUComputed_ = false;

    // Build effective viscosity and pressure gradient source
    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        if (turbulence_)
        {
            const ScalarField& nut = turbulence_->turbulentViscosity();
            nuEff_[cellIdx] = nu_ + nut[cellIdx];
        }
        else
        {
            nuEff_[cellIdx] = nu_;
        }

        // Calculate pressure gradients source term
        const Scalar volume = mesh_.cells()[cellIdx].volume();
        UxSource_[cellIdx] = -gradP_[cellIdx].x() * volume;
        UySource_[cellIdx] = -gradP_[cellIdx].y() * volume;
        UzSource_[cellIdx] = -gradP_[cellIdx].z() * volume;
    }

    // Build face-based effective viscosity
    #pragma omp parallel for schedule(static)
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            const size_t ownerIdx = face.ownerCell();

            // Check for nutWallFunction BC on wall faces
            if (turbulence_)
            {
                const BoundaryPatch& patch = face.patch()->get();
                const BoundaryData& bc =
                    bcManager_.fieldBC(patch.patchName(), Field::nut);

                if (bc.type() == BCType::nutWallFunction)
                {
                    // Use wall-function nut instead of cell-center nut
                    nuEffFace_[faceIdx] =
                        nu_ + turbulence_->nutWall()[faceIdx];
                    continue;
                }
            }

            // Other boundary faces: use owner cell value
            nuEffFace_[faceIdx] = nuEff_[ownerIdx];
        }
        else
        {
            // Internal faces: linear interpolation
            nuEffFace_[faceIdx] = interpolateToFace(face, nuEff_);
        }
    }

    // Reset interpolated diagonals accumulator
    DUf_.setAll(S(0.0));

    updateVelocityGradients();

    // Add transpose gradient source term
    if (turbulence_)
    {
        addTransposeGradientSource();
    }

    // Solve momentum equations for each component
    TransportEquation equationUx
    {
        .field          = Field::Ux,
        .phi            = Ux_,
        .flowRate       = std::cref(RhieChowFlowRatePrev_),
        .convScheme     = std::cref(momentumConvectionScheme_),
        .Gamma          = std::nullopt,
        .GammaFace      = std::cref(nuEffFace_),
        .source         = UxSource_,
        .gradPhi        = gradUx_,
        .gradScheme     = gradientScheme_
    };
    solveMomentumComponent(equationUx, UxPrev_);

    TransportEquation equationUy
    {
        .field          = Field::Uy,
        .phi            = Uy_,
        .flowRate       = std::cref(RhieChowFlowRatePrev_),
        .convScheme     = std::cref(momentumConvectionScheme_),
        .Gamma          = std::nullopt,
        .GammaFace      = std::cref(nuEffFace_),
        .source         = UySource_,
        .gradPhi        = gradUy_,
        .gradScheme     = gradientScheme_
    };
    solveMomentumComponent(equationUy, UyPrev_);

    TransportEquation equationUz
    {
        .field          = Field::Uz,
        .phi            = Uz_,
        .flowRate       = std::cref(RhieChowFlowRatePrev_),
        .convScheme     = std::cref(momentumConvectionScheme_),
        .Gamma          = std::nullopt,
        .GammaFace      = std::cref(nuEffFace_),
        .source         = UzSource_,
        .gradPhi        = gradUz_,
        .gradScheme     = gradientScheme_
    };
    solveMomentumComponent(equationUz, UzPrev_);

    // Calculate DUf_ field using complete DU_
    #pragma omp parallel for schedule(static)
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            const BoundaryData& bc =
                bcManager_.fieldBC(face.patch()->get().patchName(), Field::p);

            if (bc.type() == BCType::fixedValue)
            {
                // Fixed pressure boundary: normal pressure-velocity coupling
                DUf_[faceIdx] = DU_[face.ownerCell()];
            }
            else
            {
                // Zero gradient pressure boundary
                DUf_[faceIdx] = S(0.0);
            }
        }
        else
        {
            // Internal faces
            DUf_[faceIdx] = interpolateToFace(face, DU_);
        }
    }
}


void SIMPLE::updateRhieChowFlowRate()
{
    const size_t numFaces = mesh_.numFaces();

    #pragma omp parallel for schedule(static)
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            UxAvgf_[faceIdx] = bcManager_.boundaryFaceValue(Field::Ux, Ux_, face);
            UyAvgf_[faceIdx] = bcManager_.boundaryFaceValue(Field::Uy, Uy_, face);
            UzAvgf_[faceIdx] = bcManager_.boundaryFaceValue(Field::Uz, Uz_, face);

            const Vector Uf
            (
                UxAvgf_[faceIdx],
                UyAvgf_[faceIdx],
                UzAvgf_[faceIdx]
            );

            RhieChowFlowRate_[faceIdx] =
                dot(Uf, face.normal() * face.projectedArea());

            continue;
        }

        const size_t P = face.ownerCell();
        const size_t N = face.neighborCell().value();

        // Linear-interpolated velocity at face
        const Vector UfLinear
        (
            interpolateToFace(face, Ux_),
            interpolateToFace(face, Uy_),
            interpolateToFace(face, Uz_)
        );

        const Vector gradPAvgf = interpolateToFace(face, gradP_);
        const Vector Sf = face.normal() * face.projectedArea();
        const Vector gradPf =
            gradientScheme_.faceGradient
            (
                Field::p,
                p_,
                gradP_[P],
                gradP_[N],
                faceIdx
            );
        const Vector UfPrev
        (
            UxAvgPrevf_[faceIdx],
            UyAvgPrevf_[faceIdx],
            UzAvgPrevf_[faceIdx]
        );

        RhieChowFlowRate_[faceIdx] =
            dot(UfLinear, Sf)
          - dot((DUf_[faceIdx] * (gradPf - gradPAvgf)), Sf)
          + (S(1.0) - alphaU_)
          * (RhieChowFlowRatePrev_[faceIdx] - dot(UfPrev, Sf));
    }
}


void SIMPLE::solvePressureCorrection()
{
    const size_t numCells = mesh_.numCells();

    gradientScheme_.fieldGradient(Field::pCorr, pCorr_, gradPCorrPrecomputed_);

    // Compute mass imbalance source term
    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar net = S(0.0);
        const auto& faceIndices = mesh_.cells()[cellIdx].faceIndices();
        const auto& signs = mesh_.cells()[cellIdx].faceSigns();

        for (size_t j = 0; j < faceIndices.size(); ++j)
        {
            net += signs[j] * RhieChowFlowRate_[faceIndices[j]];
        }

        massImbalance_[cellIdx] = -net;
    }

    TransportEquation equationPCorr
    {
        .field      = Field::pCorr,
        .phi        = pCorr_,
        .flowRate   = std::nullopt,
        .convScheme = std::nullopt,
        .Gamma      = std::nullopt,
        .GammaFace  = std::cref(DUf_),
        .source     = massImbalance_,
        .gradPhi    = gradPCorrPrecomputed_,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_.buildMatrix(equationPCorr);

    // Get references to the assembled matrix
    auto& matrixA = matrixConstruct_.matrixA();
    auto& vectorB = matrixConstruct_.vectorB();

    // Map pCorr field storage as Eigen vector (zero-copy)
    pCorr_.setAll(S(0.0));

    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    pCorrSolution(pCorr_.data(), eIdx(numCells));
    pressureSolver_.solve(pCorrSolution, matrixA, vectorB);

    if (debug_)
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

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradPCorr_[cellIdx] =
            gradientScheme_.cellGradient(Field::pCorr, pCorr_, cellIdx);
    }
}


void SIMPLE::correctVelocity()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Ux_[cellIdx] -= DU_[cellIdx] * gradPCorr_[cellIdx].x();
        Uy_[cellIdx] -= DU_[cellIdx] * gradPCorr_[cellIdx].y();
        Uz_[cellIdx] -= DU_[cellIdx] * gradPCorr_[cellIdx].z();
    }

    const auto velocityConstraints =
        constraintSystem_.applyVelocityConstraints();

    if (debug_ && velocityConstraints > 0)
    {
        std::cout
            << "  Applied velocity constraints to " << velocityConstraints
            << " cells" << '\n';
    }

    // Update face velocities
    const size_t numFaces = mesh_.numFaces();

    #pragma omp parallel for schedule(static)
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            UxAvgf_[faceIdx] = bcManager_.boundaryFaceValue(Field::Ux, Ux_, face);
            UyAvgf_[faceIdx] = bcManager_.boundaryFaceValue(Field::Uy, Uy_, face);
            UzAvgf_[faceIdx] = bcManager_.boundaryFaceValue(Field::Uz, Uz_, face);
        }
        else
        {
            UxAvgf_[faceIdx] = interpolateToFace(face, Ux_);
            UyAvgf_[faceIdx] = interpolateToFace(face, Uy_);
            UzAvgf_[faceIdx] = interpolateToFace(face, Uz_);
        }
    }
}


void SIMPLE::correctFlowRate()
{
    // Update mass flux on faces
    const size_t numFaces = mesh_.numFaces();

    #pragma omp parallel for schedule(static)
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            const BoundaryData& bc =
                bcManager_.fieldBC(face.patch()->get().patchName(), Field::p);

            if
            (
                bc.type() == BCType::fixedValue
             || bc.type() == BCType::zeroGradient
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

        const size_t ownerIdx = face.ownerCell();
        const size_t neighborIdx = face.neighborCell().value();

        const Vector gradPCorrf =
            gradientScheme_.faceGradient
            (
                Field::pCorr,
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


void SIMPLE::correctPressure()
{
    Scalar sumSq = S(0.0);

    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static) reduction(+:sumSq)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        sumSq += pCorr_[cellIdx] * pCorr_[cellIdx];
    }

    lastPressureCorrectionRMS_ = std::sqrt(sumSq / S(numCells));

    // Apply pressure correction
    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        p_[cellIdx] += alphaP_ * pCorr_[cellIdx];
    }

    // Apply pressure bounds constraints
    const auto pressureConstraints =
        constraintSystem_.applyPressureConstraints();

    if (debug_ && pressureConstraints > 0)
    {
        std::cout
            << "  Applied pressure constraints to " << pressureConstraints
            << " cells" << '\n';
    }

    // Reset pressure correction for next iteration
    pCorr_.setAll(S(0.0));
}


void SIMPLE::solveTurbulence()
{
    if (turbulence_)
    {
        updateVelocityGradients();

        turbulence_->solve
        (
            Ux_,
            Uy_,
            Uz_,
            RhieChowFlowRate_,
            gradU_
        );

        lastKResidual_ = turbulence_->lastKResidual();
        lastOmegaResidual_ = turbulence_->lastOmegaResidual();
    }
}


bool SIMPLE::checkConvergence()
{
    // Compute raw residuals
    const Scalar massImbalance = this->massImbalance();
    const Scalar velocityResidual = this->velocityResidual();
    const Scalar pressureResidual = this->pressureResidual();

    // Store first-iteration references for scaling
    if (massImbalance0_ < vSmallValue)
    {
        massImbalance0_ = massImbalance;
        velocityResidual0_ = velocityResidual;
        pressureResidual0_ = pressureResidual;
        if (turbulence_)
        {
            kResidual0_ = lastKResidual_;
            omegaResidual0_ = lastOmegaResidual_;
        }
    }

    // Scale by first-iteration values
    const Scalar scaledMass = massImbalance / (massImbalance0_ + vSmallValue);

    const Scalar scaledVelocity =
        velocityResidual / (velocityResidual0_ + vSmallValue);

    const Scalar scaledPressure =
        pressureResidual / (pressureResidual0_ + vSmallValue);

    bool converged =
        (scaledMass < tolerance_)
     && (scaledVelocity < tolerance_)
     && (scaledPressure < tolerance_);

    Scalar scaledK = S(0.0);
    Scalar scaledOmega = S(0.0);

    if (turbulence_)
    {
        scaledK = lastKResidual_ / (kResidual0_ + vSmallValue);
        scaledOmega = lastOmegaResidual_ / (omegaResidual0_ + vSmallValue);

        converged = converged
            && (scaledK < tolerance_)
            && (scaledOmega < tolerance_);
    }

    if (debug_)
    {
        Logger::subsection("Scaled residuals");
        Logger::scaledResidual("mass",     scaledMass);
        Logger::scaledResidual("velocity", scaledVelocity);
        Logger::scaledResidual("pressure", scaledPressure);
        if (turbulence_)
        {
            Logger::scaledResidual("k",     scaledK);
            Logger::scaledResidual("omega", scaledOmega);
        }
    }
    else if (turbulence_)
    {
        Logger::residualSummary
        (
            scaledMass,
            scaledVelocity,
            scaledPressure,
            scaledK,
            scaledOmega
        );
    }
    else
    {
        Logger::residualSummary(scaledMass, scaledVelocity, scaledPressure);
    }

    return converged;
}


void SIMPLE::updateVelocityGradients()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
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

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
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


void SIMPLE::addTransposeGradientSource()
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar sumX = S(0.0);
        Scalar sumY = S(0.0);
        Scalar sumZ = S(0.0);

        const auto& cell = mesh_.cells()[cellIdx];
        const auto faceIndices = cell.faceIndices();
        const auto faceSigns = cell.faceSigns();

        for (size_t j = 0; j < faceIndices.size(); ++j)
        {
            const size_t faceIdx = faceIndices[j];
            const Scalar sign = S(faceSigns[j]);
            const Face& face = mesh_.faces()[faceIdx];

            const Vector Sf =
                face.normal() * face.projectedArea() * sign;
            const Scalar nuEfff = nuEffFace_[faceIdx];

            Tensor gradUf;
            if (face.isBoundary())
            {
                gradUf = gradU_[cellIdx];
            }
            else
            {
                gradUf = interpolateToFace(face, gradU_);
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


void SIMPLE::solveMomentumComponent
(
    TransportEquation& eq,
    const ScalarField& componentPrev
)
{
    matrixConstruct_.buildMatrix(eq);

    auto& matrixA = matrixConstruct_.matrixA();
    auto& vectorB = matrixConstruct_.vectorB();

    matrixConstruct_.relax(alphaU_, componentPrev);

    const size_t numCells = mesh_.numCells();

    // DU_ is identical for all 3 momentum components — compute only once
    if (DUComputed_ == false)
    {
        #pragma omp parallel for schedule(static)
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            DU_[cellIdx] =
                mesh_.cells()[cellIdx].volume()
              / (matrixA.coeff(eIdx(cellIdx), eIdx(cellIdx)) + vSmallValue);
        }
        DUComputed_ = true;
    }

    // Map phi directly as Eigen vector: the zero-copy solve writes the
    // result straight into the bound velocity component field
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    solutionMap(eq.phi.data(), eIdx(numCells));

    momentumSolver_.solve(solutionMap, matrixA, vectorB);

    if (debug_)
    {
        const SolvePerformance& momentumPerformance =
            momentumSolver_.lastPerformance();

        Logger::residualRow
        (
            fieldToString(eq.field),
            momentumPerformance.solverName,
            momentumPerformance.iterations,
            momentumPerformance.finalResidual
        );
    }
}


Scalar SIMPLE::massImbalance() const noexcept
{
    // Dimensionless normalized continuity residual per cell, averaged
    Scalar totalNormImbalance = S(0.0);

    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static) reduction(+:totalNormImbalance)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar net = S(0.0);
        Scalar sumAbs = S(0.0);

        for (size_t j = 0; j < mesh_.cells()[cellIdx].faceIndices().size(); ++j)
        {
            const size_t faceIdx = mesh_.cells()[cellIdx].faceIndices()[j];
            const int sign = mesh_.cells()[cellIdx].faceSigns()[j];
            const Scalar mf = RhieChowFlowRate_[faceIdx];
            net += S(sign) * mf;
            sumAbs += std::abs(mf);
        }

        const Scalar denom = sumAbs + vSmallValue;
        totalNormImbalance += std::abs(net) / denom;
    }

    return totalNormImbalance / S(std::max<size_t>(1, numCells));
}


Scalar SIMPLE::velocityResidual() const noexcept
{
    // Normalized residual: ||U - U_prev||_2 / (||U_prev||_2 + eps)
    Scalar num = S(0.0);
    Scalar den = S(0.0);

    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static) reduction(+:num, den)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar dx = Ux_[cellIdx] - UxPrev_[cellIdx];
        const Scalar dy = Uy_[cellIdx] - UyPrev_[cellIdx];
        const Scalar dz = Uz_[cellIdx] - UzPrev_[cellIdx];

        num += dx * dx + dy * dy + dz * dz;
        den += UxPrev_[cellIdx] * UxPrev_[cellIdx]
             + UyPrev_[cellIdx] * UyPrev_[cellIdx]
             + UzPrev_[cellIdx] * UzPrev_[cellIdx];
    }

    num = std::sqrt(num + vSmallValue);
    den = std::sqrt(den + vSmallValue);

    return num / den;
}


Scalar SIMPLE::pressureResidual() const noexcept
{
    // Normalize p' RMS by RMS(p)
    Scalar sumP2 = S(0.0);

    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static) reduction(+:sumP2)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        sumP2 += p_[cellIdx] * p_[cellIdx];
    }

    const Scalar pRms = std::sqrt(sumP2 / S(numCells));

    return lastPressureCorrectionRMS_ / (pRms + vSmallValue);
}
