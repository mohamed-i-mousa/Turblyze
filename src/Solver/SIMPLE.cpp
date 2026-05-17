/******************************************************************************
 * @file SIMPLE.cpp
 * @brief Implementation of the SIMPLE algorithm for pressure-velocity coupling
 *****************************************************************************/

#include "SIMPLE.h"

#include <cmath>
#include <iostream>
#include <algorithm>

#include <omp.h>

#include "LinearInterpolation.h"
#include "Logger.h"
#include "Scalar.h"
#include "kOmegaSST.h"
#include "Constraint.h"

// ************************* Special Member Functions *************************

SIMPLE::SIMPLE
(
    const Mesh& mesh,
    const BoundaryConditions& bc,
    const GradientScheme& gradScheme,
    const ConvectionScheme& convSchemes
) 
:
    mesh_(mesh),
    bcManager_(bc),
    gradientScheme_(gradScheme),
    convectionScheme_(convSchemes)
{}

// ****************************** Setter Methods ******************************

void SIMPLE::setTurbulenceSolvers
(
    const LinearSolver& kSolver,
    const LinearSolver& omegaSolver
)
{
    if (turbulenceModel_)
    {
        turbulenceModel_->kSolverSettings() = kSolver;
        turbulenceModel_->omegaSolverSettings() = omegaSolver;
    }
}

// ******************************* Public Methods ******************************

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


void SIMPLE::initialize
(
    const Vector& initialVelocity,
    Scalar initialPressure,
    Scalar initialK,
    Scalar initialOmega,
    bool enableTurbulence
)
{
    matrixConstruct_ = std::make_unique<Matrix>(mesh_, bcManager_);
    constraintSystem_ = std::make_unique<Constraint>(Ux_, Uy_, Uz_, p_);

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

    if (enableTurbulence)
    {
        turbulenceModel_ =
            std::make_unique<kOmegaSST>
            (
                mesh_,
                bcManager_,
                gradientScheme_,
                convectionScheme_.k(),
                convectionScheme_.omega()
            );

        turbulenceModel_->setDebug(debug_);

        turbulenceModel_->
            initialize
            (
                nu_,
                initialK,
                initialOmega,
                alphaK_,
                alphaOmega_
            );
    }
}


void SIMPLE::updateVelocityGradients
(
    VectorField& gradUx,
    VectorField& gradUy,
    VectorField& gradUz
)
{
    const size_t numCells = mesh_.numCells();

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradUx[cellIdx] =
            gradientScheme_.cellGradient(Field::Ux, Ux_, cellIdx);
        gradUy[cellIdx] =
            gradientScheme_.cellGradient(Field::Uy, Uy_, cellIdx);
        gradUz[cellIdx] =
            gradientScheme_.cellGradient(Field::Uz, Uz_, cellIdx);
    }

    gradientScheme_.limitGradient(Field::Ux, Ux_, gradUx);
    gradientScheme_.limitGradient(Field::Uy, Uy_, gradUy);
    gradientScheme_.limitGradient(Field::Uz, Uz_, gradUz);

    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradU_[cellIdx] =
            Tensor::fromRows
            (
                gradUx[cellIdx],
                gradUy[cellIdx],
                gradUz[cellIdx]
            );
    }
}


void SIMPLE::solveMomentumEquations()
{
    const size_t numCells = mesh_.numCells();
    const size_t numFaces = mesh_.numFaces();

    ScalarField nuEff;
    ScalarField UxSource;
    ScalarField UySource;
    ScalarField UzSource;

    // Reset diagonals accumulator
    DU_.setAll(S(0.0));
    DUComputed_ = false;

    // Build effective viscosity and pressure gradient source
    #pragma omp parallel for schedule(static)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        if (turbulenceModel_)
        {
            const ScalarField& nut = turbulenceModel_->turbulentViscosity();
            nuEff[cellIdx] = nu_ + nut[cellIdx];
        }
        else
        {
            nuEff[cellIdx] = nu_;
        }

        // Calculate pressure gradients source term
        UxSource[cellIdx] = -gradP_[cellIdx].x() * mesh_.cells()[cellIdx].volume();
        UySource[cellIdx] = -gradP_[cellIdx].y() * mesh_.cells()[cellIdx].volume();
        UzSource[cellIdx] = -gradP_[cellIdx].z() * mesh_.cells()[cellIdx].volume();
    }

    // Build face-based effective viscosity
    FaceData<Scalar> nuEffFace(nu_);

    #pragma omp parallel for schedule(static)
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            const size_t ownerIdx = face.ownerCell();

            // Check for NUT_WALL_FUNCTION BC on wall faces
            if (turbulenceModel_)
            {
                const BoundaryPatch& patch = face.patch()->get();
                const BoundaryData& bc =
                    bcManager_.fieldBC(patch.patchName(), Field::nut);

                if (bc.type() == BCType::NUT_WALL_FUNCTION)
                {
                    // Use wall-function nut instead of cell-center nut
                    nuEffFace[faceIdx] =
                        nu_ + turbulenceModel_->nutWall()[faceIdx];
                    continue;
                }
            }

            // Other boundary faces: use owner cell value
            nuEffFace[faceIdx] = nuEff[ownerIdx];
        }
        else
        {
            // Internal faces: linear interpolation
            nuEffFace[faceIdx] = interpolateToFace(face, nuEff);
        }
    }

    // Reset interpolated diagonals accumulator
    DUf_.setAll(S(0.0));

    VectorField gradUx;
    VectorField gradUy;
    VectorField gradUz;
    updateVelocityGradients(gradUx, gradUy, gradUz);

    // Add transpose gradient source term
    if (turbulenceModel_)
    {
        ScalarField transposeSourceX;
        ScalarField transposeSourceY;
        ScalarField transposeSourceZ;

        transposeGradientSource
        (
            nuEffFace,
            transposeSourceX,
            transposeSourceY,
            transposeSourceZ
        );

        // Add transpose gradient contribution to momentum source terms
        #pragma omp parallel for schedule(static)
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            UxSource[cellIdx] += transposeSourceX[cellIdx];
            UySource[cellIdx] += transposeSourceY[cellIdx];
            UzSource[cellIdx] += transposeSourceZ[cellIdx];
        }
    }

    // Solve momentum equations for each component
    TransportEquation equationUx
    {
        .field          = Field::Ux,
        .phi            = Ux_,
        .flowRate       = std::cref(RhieChowFlowRatePrev_),
        .convScheme     = std::cref(convectionScheme_.momentum()),
        .Gamma          = std::nullopt,
        .GammaFace      = std::cref(nuEffFace),
        .source         = UxSource,
        .gradPhi        = gradUx,
        .gradScheme     = gradientScheme_
    };
    solveMomentumEquation(equationUx, UxPrev_);

    TransportEquation equationUy
    {
        .field          = Field::Uy,
        .phi            = Uy_,
        .flowRate       = std::cref(RhieChowFlowRatePrev_),
        .convScheme     = std::cref(convectionScheme_.momentum()),
        .Gamma          = std::nullopt,
        .GammaFace      = std::cref(nuEffFace),
        .source         = UySource,
        .gradPhi        = gradUy,
        .gradScheme     = gradientScheme_
    };
    solveMomentumEquation(equationUy, UyPrev_);

    TransportEquation equationUz
    {
        .field          = Field::Uz,
        .phi            = Uz_,
        .flowRate       = std::cref(RhieChowFlowRatePrev_),
        .convScheme     = std::cref(convectionScheme_.momentum()),
        .Gamma          = std::nullopt,
        .GammaFace      = std::cref(nuEffFace),
        .source         = UzSource,
        .gradPhi        = gradUz,
        .gradScheme     = gradientScheme_
    };
    solveMomentumEquation(equationUz, UzPrev_);

    // Calculate DUf_ field using complete DU_
    #pragma omp parallel for schedule(static)
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            const BoundaryData& bc =
                bcManager_.fieldBC(face.patch()->get().patchName(), Field::p);

            if (bc.type() == BCType::FIXED_VALUE)
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

    VectorField gradPCorrPrecomputed;
    gradientScheme_.fieldGradient(Field::pCorr, pCorr_, gradPCorrPrecomputed);

    // Compute mass imbalance source term
    ScalarField massImbalance;

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

        massImbalance[cellIdx] = -net;
    }

    TransportEquation equationPCorr
    {
        .field      = Field::pCorr,
        .phi        = pCorr_,
        .flowRate   = std::nullopt,
        .convScheme = std::nullopt,
        .Gamma      = std::nullopt,
        .GammaFace  = std::cref(DUf_),
        .source     = massImbalance,
        .gradPhi    = gradPCorrPrecomputed,
        .gradScheme = gradientScheme_
    };

    matrixConstruct_->buildMatrix(equationPCorr);

    // Get references to the assembled matrix
    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    // Map pCorr field storage as Eigen vector (zero-copy)
    pCorr_.setAll(S(0.0));

    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    pCorrSolution(pCorr_.data(), eIdx(numCells));
    pressureSolver_.solveWithPCG(pCorrSolution, matrixA, vectorB);

    if (debug_)
    {
        Logger::residualRow
        (
            "p'",
            "PCG",
            pressureSolver_.lastIterations(),
            pressureSolver_.lastResidual()
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
        constraintSystem_->applyVelocityConstraints();

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
        constraintSystem_->applyPressureConstraints();

    if (debug_ && pressureConstraints > 0)
    {
        std::cout
            << "  Applied pressure constraints to " << pressureConstraints
            << " cells" << '\n';
    }

    // Reset pressure correction for next iteration
    pCorr_.setAll(S(0.0));
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
                bc.type() == BCType::FIXED_VALUE
             || bc.type() == BCType::ZERO_GRADIENT
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


void SIMPLE::solveTurbulence()
{
    if (turbulenceModel_)
    {
        const size_t numCells = mesh_.numCells();

        VectorField gradUx;
        VectorField gradUy;
        VectorField gradUz;
        updateVelocityGradients(gradUx, gradUy, gradUz);

        const auto& kField = turbulenceModel_->k();
        const auto& omegaField = turbulenceModel_->omega();

        ScalarField kPrev = kField;
        ScalarField omegaPrev = omegaField;

        turbulenceModel_->solve(Ux_, Uy_, Uz_, RhieChowFlowRate_, gradU_);

        // Compute normalised change: ||x - x_prev|| / ||x_prev||
        Scalar kDiffSq = S(0.0);
        Scalar kPrevSq = S(0.0);
        Scalar omDiffSq = S(0.0);
        Scalar omPrevSq = S(0.0);

        #pragma omp parallel for schedule(static) \
            reduction(+:kDiffSq, kPrevSq, omDiffSq, omPrevSq)
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            const Scalar dk = kField[cellIdx] - kPrev[cellIdx];
            kDiffSq += dk * dk;
            kPrevSq += kPrev[cellIdx] * kPrev[cellIdx];

            const Scalar dw = omegaField[cellIdx] - omegaPrev[cellIdx];
            omDiffSq += dw * dw;
            omPrevSq += omegaPrev[cellIdx] * omegaPrev[cellIdx];
        }

        lastKResidual_ =
            std::sqrt(kDiffSq + vSmallValue)
          / std::sqrt(kPrevSq + vSmallValue);
        lastOmegaResidual_ =
            std::sqrt(omDiffSq + vSmallValue)
          / std::sqrt(omPrevSq + vSmallValue);
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
        if (turbulenceModel_)
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

    if (turbulenceModel_)
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
        if (turbulenceModel_)
        {
            Logger::scaledResidual("k",     scaledK);
            Logger::scaledResidual("omega", scaledOmega);
        }
    }
    else if (turbulenceModel_)
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

// ****************************** Private Methods ******************************

Scalar SIMPLE::massImbalance() const
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


Scalar SIMPLE::velocityResidual() const
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


Scalar SIMPLE::pressureResidual() const
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


void SIMPLE::transposeGradientSource
(
    const FaceData<Scalar>& nuEffFace,
    ScalarField& transposeSourceX,
    ScalarField& transposeSourceY,
    ScalarField& transposeSourceZ
) const
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
            const Scalar nuEfff = nuEffFace[faceIdx];

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

        transposeSourceX[cellIdx] = sumX;
        transposeSourceY[cellIdx] = sumY;
        transposeSourceZ[cellIdx] = sumZ;
    }
}


void SIMPLE::solveMomentumEquation
(
    TransportEquation& eq,
    const ScalarField& componentPrev
)
{
    matrixConstruct_->buildMatrix(eq);

    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    matrixConstruct_->relax(alphaU_, componentPrev);

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

    momentumSolver_.solveWithBiCGSTAB
    (
        solutionMap,
        matrixA,
        vectorB
    );

    if (debug_)
    {
        Logger::residualRow
        (
            fieldToString(eq.field),
            "BiCGSTAB",
            momentumSolver_.lastIterations(),
            momentumSolver_.lastResidual()
        );
    }
}
