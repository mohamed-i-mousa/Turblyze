/******************************************************************************
 * @file SIMPLE.cpp
 * @brief Implementation of the SIMPLE algorithm for pressure-velocity coupling
 *****************************************************************************/

#include "SIMPLE.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

#include "LinearInterpolation.hpp"
#include "Scalar.hpp"
#include "kOmegaSST.hpp"
#include "Constraint.hpp"


// ************************ Constructor & Destructor ************************

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
    std::cout
        << "\n=== Starting SIMPLE Loop ===" << std::endl
        << std::endl;

    // Reset first-iteration residual references for clean convergence tracking
    massImbalance0_    = S(0.0);
    velocityResidual0_ = S(0.0);
    pressureResidual0_ = S(0.0);
    kResidual0_        = S(0.0);
    omegaResidual0_    = S(0.0);

    int iteration = 0;

    bool converged = false;

    while (!converged && iteration < maxIterations_)
    {
        std::cout
            << " Iteration " << iteration + 1;

        UPrev_ = U_;
        UAvgPrevf_ = UAvgf_;
        RhieChowFlowRatePrev_ = RhieChowFlowRate_;

        size_t numCells = mesh_.numCells();
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            gradP_[cellIdx] = gradientScheme_.cellGradient("p", p_, cellIdx);
        }

        solveMomentumEquations();

        updateRhieChowFlowRate();

        solvePressureCorrection();

        correctVelocity();

        correctFlowRate();

        correctPressure();

        solveTurbulence();

        converged = checkConvergence();

        iteration++;
    }

    if (!converged)
    {
        std::cout
            << "WARNING: SIMPLE algorithm did not converge after "
            << maxIterations_ << " iterations." << std::endl;
    }
    else
    {
        std::cout
            << "SIMPLE algorithm converged in " << iteration
            << " iterations." << std::endl;
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

    // Initialize constraint system
    constraintSystem_ = std::make_unique<Constraint>(U_, p_);

    U_.setAll(initialVelocity);
    p_.setAll(initialPressure);
    UAvgf_.setAll(initialVelocity);

    // Initialize RhieChowFlowRate_ with linear interpolation
    size_t numFaces = mesh_.numFaces();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        Vector Uf;
        if (face.isBoundary())
        {
            Uf = bcManager_.boundaryVectorFaceValue("U", U_, face);
        }
        else
        {
            Uf = interpolateToFace(face, U_);
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

        if (debug_)
        {
            std::cout
                << "k-omega SST turbulence model initialized." << std::endl;
        }
    }

    if (debug_)
    {
        std::cout
            << "SIMPLE algorithm is initialized!" << std::endl;
    }
}

void SIMPLE::solveMomentumEquations()
{
    size_t numCells = mesh_.numCells();
    size_t numFaces = mesh_.numFaces();

    ScalarField nuEff;

    // Extract velocity components
    ScalarField Ux = extractComponent(U_, Component::X);
    ScalarField Uy = extractComponent(U_, Component::Y);
    ScalarField Uz = extractComponent(U_, Component::Z);

    // Extract previous-iteration components
    ScalarField UxPrev = extractComponent(UPrev_, Component::X);
    ScalarField UyPrev = extractComponent(UPrev_, Component::Y);
    ScalarField UzPrev = extractComponent(UPrev_, Component::Z);

    ScalarField UxSource;
    ScalarField UySource;
    ScalarField UzSource;

    // Reset diagonals accumulator
    DU_.setAll(0.0);

    // Build effective viscosity and pressure gradient source
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

    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            size_t ownerIdx = face.ownerCell();

            // Check for NUT_WALL_FUNCTION BC on wall faces
            if (turbulenceModel_)
            {
                const BoundaryPatch& patch = face.patch()->get();
                const BoundaryData& bc =
                    bcManager_.fieldBC(patch.patchName(), "nut");

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
    DUf_.setAll(0.0);

    // Compute per-row velocity gradients and assemble the tensor field.
    // Row VectorFields feed the scalar momentum TransportEquations;
    // gradU_ (TensorField) is reused by the transpose source and the
    // turbulence model.
    VectorField gradUx_row;
    VectorField gradUy_row;
    VectorField gradUz_row;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradUx_row[cellIdx] =
            gradientScheme_.cellGradient("U", Ux, cellIdx, nullptr, 0);
        gradUy_row[cellIdx] =
            gradientScheme_.cellGradient("U", Uy, cellIdx, nullptr, 1);
        gradUz_row[cellIdx] =
            gradientScheme_.cellGradient("U", Uz, cellIdx, nullptr, 2);

        gradU_[cellIdx] = Tensor::fromRows
        (
            gradUx_row[cellIdx],
            gradUy_row[cellIdx],
            gradUz_row[cellIdx]
        );
    }

    // Add transpose gradient source term
    if (turbulenceModel_)
    {
        ScalarField transposeSourceX;
        ScalarField transposeSourceY;
        ScalarField transposeSourceZ;

        calculateTransposeGradientSource
        (
            nuEffFace,
            transposeSourceX,
            transposeSourceY,
            transposeSourceZ
        );

        // Add transpose gradient contribution to momentum source terms
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
        .fieldName      = "U",
        .phi            = Ux,
        .flowRate       = std::cref(RhieChowFlowRatePrev_),
        .convScheme     = std::cref(convectionScheme_.momentum()),
        .Gamma          = std::nullopt,
        .GammaFace      = std::cref(nuEffFace),
        .source         = UxSource,
        .gradPhi        = gradUx_row,
        .gradScheme     = gradientScheme_,
        .componentIdx   = 0
    };
    solveMomentumComponent(Component::X, equationUx, UxPrev);

    TransportEquation equationUy
    {
        .fieldName      = "U",
        .phi            = Uy,
        .flowRate       = std::cref(RhieChowFlowRatePrev_),
        .convScheme     = std::cref(convectionScheme_.momentum()),
        .Gamma          = std::nullopt,
        .GammaFace      = std::cref(nuEffFace),
        .source         = UySource,
        .gradPhi        = gradUy_row,
        .gradScheme     = gradientScheme_,
        .componentIdx   = 1
    };
    solveMomentumComponent(Component::Y, equationUy, UyPrev);

    TransportEquation equationUz
    {
        .fieldName      = "U",
        .phi            = Uz,
        .flowRate       = std::cref(RhieChowFlowRatePrev_),
        .convScheme     = std::cref(convectionScheme_.momentum()),
        .Gamma          = std::nullopt,
        .GammaFace      = std::cref(nuEffFace),
        .source         = UzSource,
        .gradPhi        = gradUz_row,
        .gradScheme     = gradientScheme_,
        .componentIdx   = 2
    };
    solveMomentumComponent(Component::Z, equationUz, UzPrev);

    // Calculate DUf_ field using complete DU_
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            const BoundaryData& bc =
                bcManager_.fieldBC(face.patch()->get().patchName(), "p");

            if (bc.type() == BCType::FIXED_VALUE)
            {
                // Fixed pressure boundary: normal pressure-velocity coupling
                DUf_[faceIdx] = DU_[face.ownerCell()];
            }
            else
            {
                // Zero gradient pressure boundary
                DUf_[faceIdx] = 0.0;
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
    size_t numFaces = mesh_.numFaces();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            UAvgf_[faceIdx] =
                bcManager_.boundaryVectorFaceValue("U", U_, face);

            RhieChowFlowRate_[faceIdx] =
                dot(UAvgf_[faceIdx], face.normal() * face.projectedArea());

            continue;
        }

        const size_t P = face.ownerCell();
        const size_t N = face.neighborCell().value();

        // Linear-interpolated velocity at face
        const Vector UfLinear = interpolateToFace(face, U_);

        const Vector gradPAvgf = interpolateToFace(face, gradP_);

        const Vector Sf = face.normal() * face.projectedArea();

        Vector gradPf =
            gradientScheme_.faceGradient
            (
                "p",
                p_,
                gradP_[P],
                gradP_[N],
                faceIdx
            );

        RhieChowFlowRate_[faceIdx] =
            dot(UfLinear, Sf)
          - dot((DUf_[faceIdx] * (gradPf - gradPAvgf)), Sf)
          + (S(1.0) - alphaU_)
          * (RhieChowFlowRatePrev_[faceIdx] - dot(UAvgPrevf_[faceIdx], Sf));
    }
}

void SIMPLE::solvePressureCorrection()
{
    size_t numCells = mesh_.numCells();

    VectorField gradPCorrPrecomputed;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradPCorrPrecomputed[cellIdx] =
            gradientScheme_.cellGradient("pCorr", pCorr_, cellIdx);
    }

    // Compute mass imbalance source term
    ScalarField massImbalance;

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
        .fieldName  = "pCorr",
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
    pCorrSolution(pCorr_.data(), static_cast<Eigen::Index>(numCells));

    pressureSolver_.solveWithPCG(pCorrSolution, matrixA, vectorB);

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradPCorr_[cellIdx] =
            gradientScheme_.cellGradient("pCorr", pCorr_, cellIdx);
    }
}

void SIMPLE::correctVelocity()
{
    size_t numCells = mesh_.numCells();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        U_[cellIdx].setX
        (
            U_[cellIdx].x() - DU_[cellIdx] * gradPCorr_[cellIdx].x()
        );

        U_[cellIdx].setY
        (
            U_[cellIdx].y() - DU_[cellIdx] * gradPCorr_[cellIdx].y()
        );

        U_[cellIdx].setZ
        (
            U_[cellIdx].z() - DU_[cellIdx] * gradPCorr_[cellIdx].z()
        );
    }

    auto velocityConstraints =
        constraintSystem_->applyVelocityConstraints();

    if (debug_ && velocityConstraints > 0)
    {
        std::cout
            << "  Applied velocity constraints to " << velocityConstraints
            << " cells" << std::endl;
    }

    // Update face velocities
    size_t numFaces = mesh_.numFaces();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            UAvgf_[faceIdx] =
                bcManager_.boundaryVectorFaceValue("U", U_, face);
        }
        else
        {
            UAvgf_[faceIdx] = interpolateToFace(face, U_);
        }
    }
}

void SIMPLE::correctPressure()
{
    Scalar sumSq = 0.0;

    size_t numCells = mesh_.numCells();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        sumSq += pCorr_[cellIdx] * pCorr_[cellIdx];
    }

    lastPressureCorrectionRMS_ = std::sqrt(sumSq / S(numCells));

    // Apply pressure correction
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        p_[cellIdx] += alphaP_ * pCorr_[cellIdx];
    }

    // Apply pressure bounds constraints
    auto pressureConstraints =
        constraintSystem_->applyPressureConstraints();

    if (debug_ && pressureConstraints > 0)
    {
        std::cout
            << "  Applied pressure constraints to " << pressureConstraints
            << " cells" << std::endl;
    }

    // Reset pressure correction for next iteration
    pCorr_.setAll(0.0);
}

void SIMPLE::correctFlowRate()
{
    // Update mass flux on faces
    size_t numFaces = mesh_.numFaces();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

        if (face.isBoundary())
        {
            const BoundaryData& bc =
                bcManager_.fieldBC(face.patch()->get().patchName(), "p");

            if
            (
                bc.type() == BCType::FIXED_VALUE
             || bc.type() == BCType::ZERO_GRADIENT
            )
            {
                continue;
            }

            Scalar gradn = dot(gradPCorr_[face.ownerCell()], face.normal());

            Scalar flowRateCorrection =
                DU_[face.ownerCell()] * gradn * face.projectedArea();

            RhieChowFlowRate_[faceIdx] -= flowRateCorrection;
            continue;
        }

        size_t ownerIdx = face.ownerCell();
        size_t neighborIdx = face.neighborCell().value();

        Vector gradPCorrf =
            gradientScheme_.faceGradient
            (
                "pCorr",
                pCorr_,
                gradPCorr_[ownerIdx],
                gradPCorr_[neighborIdx],
                faceIdx
            );

        Vector Sf = face.normal() * face.projectedArea();

        Scalar flowRateCorrection = DUf_[faceIdx] * dot(gradPCorrf, Sf);

        RhieChowFlowRate_[faceIdx] -= flowRateCorrection;
    }
}

void SIMPLE::solveTurbulence()
{
    if (turbulenceModel_)
    {
        if (debug_)
        {
            std::cout
                << "  Solving turbulence equations..."
                << std::endl;
        }

        // Recompute velocity gradients from corrected velocity
        size_t numCells = mesh_.numCells();

        ScalarField Ux = extractComponent(U_, Component::X);
        ScalarField Uy = extractComponent(U_, Component::Y);
        ScalarField Uz = extractComponent(U_, Component::Z);

        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            Vector gUx =
                gradientScheme_.cellGradient("U", Ux, cellIdx, nullptr, 0);
            Vector gUy =
                gradientScheme_.cellGradient("U", Uy, cellIdx, nullptr, 1);
            Vector gUz =
                gradientScheme_.cellGradient("U", Uz, cellIdx, nullptr, 2);

            gradU_[cellIdx] = Tensor::fromRows(gUx, gUy, gUz);
        }

        const auto& kField = turbulenceModel_->k();
        const auto& omegaField = turbulenceModel_->omega();

        ScalarField kPrev;
        ScalarField omegaPrev;
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            kPrev[cellIdx] = kField[cellIdx];
            omegaPrev[cellIdx] = omegaField[cellIdx];
        }

        turbulenceModel_->solve(U_, RhieChowFlowRate_, gradU_);

        // Compute normalised change: ||x - x_prev|| / ||x_prev||
        Scalar kDiffSq = 0.0;
        Scalar kPrevSq = 0.0;
        Scalar omDiffSq = 0.0;
        Scalar omPrevSq = 0.0;
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            Scalar dk = kField[cellIdx] - kPrev[cellIdx];
            kDiffSq += dk * dk;
            kPrevSq += kPrev[cellIdx] * kPrev[cellIdx];

            Scalar dw = omegaField[cellIdx] - omegaPrev[cellIdx];
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
    Scalar massImbalance = this->massImbalance();
    Scalar velocityResidual = this->velocityResidual();
    Scalar pressureResidual = this->pressureResidual();

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
    Scalar scaledMass = massImbalance
        / (massImbalance0_ + vSmallValue);
    Scalar scaledVelocity = velocityResidual
        / (velocityResidual0_ + vSmallValue);
    Scalar scaledPressure = pressureResidual
        / (pressureResidual0_ + vSmallValue);

    bool converged =
        (scaledMass < tolerance_)
     && (scaledVelocity < tolerance_)
     && (scaledPressure < tolerance_);

    if (turbulenceModel_)
    {
        Scalar scaledK = lastKResidual_
            / (kResidual0_ + vSmallValue);
        Scalar scaledOmega = lastOmegaResidual_
            / (omegaResidual0_ + vSmallValue);

        converged = converged
            && (scaledK < tolerance_)
            && (scaledOmega < tolerance_);

        std::cout
            << " - Mass: " << std::scientific
            << scaledMass
            << ", Velocity: " << scaledVelocity
            << ", Pressure: " << scaledPressure
            << ", k: " << scaledK
            << ", omega: " << scaledOmega
            << std::fixed << std::endl;
    }
    else
    {
        std::cout
            << " - Mass: " << std::scientific
            << scaledMass
            << ", Velocity: " << scaledVelocity
            << ", Pressure: " << scaledPressure
            << std::fixed << std::endl;
    }

    return converged;
}


// ****************************** Private Methods ******************************

ScalarField SIMPLE::extractComponent
(
    const VectorField& V,
    Component component
)
{
    size_t numCells = V.size();

    ScalarField result;

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        switch (component)
        {
            case Component::X: result[cellIdx] = V[cellIdx].x(); break;
            case Component::Y: result[cellIdx] = V[cellIdx].y(); break;
            case Component::Z: result[cellIdx] = V[cellIdx].z(); break;
        }
    }

    return result;
}

Scalar SIMPLE::massImbalance() const
{
    // Dimensionless normalized continuity residual per cell, averaged
    Scalar totalNormImbalance = 0.0;

    size_t numCells = mesh_.numCells();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        Scalar net = 0.0;
        Scalar sumAbs = 0.0;

        for (size_t j = 0; j < mesh_.cells()[cellIdx].faceIndices().size(); ++j)
        {
            const size_t faceIdx = mesh_.cells()[cellIdx].faceIndices()[j];
            const int sign = mesh_.cells()[cellIdx].faceSigns()[j];
            const Scalar mf = RhieChowFlowRate_[faceIdx];
            net += sign * mf;
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
    Scalar num = 0.0;
    Scalar den = 0.0;

    size_t numCells = mesh_.numCells();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Vector d = U_[cellIdx] - UPrev_[cellIdx];
        num += d.magnitudeSquared();
        den += UPrev_[cellIdx].magnitudeSquared();
    }

    num = std::sqrt(num + vSmallValue);
    den = std::sqrt(den + vSmallValue);

    return num / den;
}

Scalar SIMPLE::pressureResidual() const
{
    // Normalize p' RMS by RMS(p)
    Scalar sumP2 = 0.0;

    size_t numCells = mesh_.numCells();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        sumP2 += p_[cellIdx] * p_[cellIdx];
    }

    const Scalar pRms = std::sqrt(sumP2 / S(numCells));

    return lastPressureCorrectionRMS_ / (pRms + vSmallValue);
}

void SIMPLE::calculateTransposeGradientSource
(
    const FaceData<Scalar>& nuEffFace,
    ScalarField& transposeSourceX,
    ScalarField& transposeSourceY,
    ScalarField& transposeSourceZ
)
{
    transposeSourceX.setAll(0.0);
    transposeSourceY.setAll(0.0);
    transposeSourceZ.setAll(0.0);

    for (size_t faceIdx = 0; faceIdx < mesh_.numFaces(); ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];
        const size_t ownerIdx = face.ownerCell();
        Vector Sf = face.normal() * face.projectedArea();

        Scalar nuEfff = nuEffFace[faceIdx];

        // Columns of the transposed velocity-gradient tensor are the
        // rows of grad(U). Flux for component j is:
        //   nuEff_f * dot(gradU_f.col(j), Sf)
        if (face.isBoundary())
        {
            // For boundary faces, use the owner cell tensor directly
            const Tensor& gradUf = gradU_[ownerIdx];

            transposeSourceX[ownerIdx] += nuEfff * dot(gradUf.col(0), Sf);
            transposeSourceY[ownerIdx] += nuEfff * dot(gradUf.col(1), Sf);
            transposeSourceZ[ownerIdx] += nuEfff * dot(gradUf.col(2), Sf);
        }
        else
        {
            // Internal face: interpolate tensor to face
            const size_t neighborIdx = face.neighborCell().value();

            Tensor gradUf = interpolateToFace(face, gradU_);

            Scalar fluxX = nuEfff * dot(gradUf.col(0), Sf);
            Scalar fluxY = nuEfff * dot(gradUf.col(1), Sf);
            Scalar fluxZ = nuEfff * dot(gradUf.col(2), Sf);

            transposeSourceX[ownerIdx] += fluxX;
            transposeSourceY[ownerIdx] += fluxY;
            transposeSourceZ[ownerIdx] += fluxZ;

            transposeSourceX[neighborIdx] -= fluxX;
            transposeSourceY[neighborIdx] -= fluxY;
            transposeSourceZ[neighborIdx] -= fluxZ;
        }
    }
}

void SIMPLE::solveMomentumComponent
(
    Component component,
    TransportEquation& eq,
    const ScalarField& componentPrev
)
{
    matrixConstruct_->buildMatrix(eq);

    auto& matrixA = matrixConstruct_->matrixA();
    auto& vectorB = matrixConstruct_->vectorB();

    matrixConstruct_->relax(alphaU_, componentPrev);

    size_t numCells = mesh_.numCells();

    // DU_ is identical for all 3 momentum components — compute only once
    if (component == Component::X)
    {
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            DU_[cellIdx] =
                mesh_.cells()[cellIdx].volume()
              / (matrixA.coeff(static_cast<Eigen::Index>(cellIdx), static_cast<Eigen::Index>(cellIdx)) + vSmallValue);
        }
    }

    // Map phi directly as Eigen vector (zero-copy solve)
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    solutionMap(eq.phi.data(), static_cast<Eigen::Index>(numCells));

    momentumSolver_.solveWithBiCGSTAB
    (
        solutionMap,
        matrixA,
        vectorB
    );

    // Write solved component back to velocity VectorField
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        switch (component)
        {
            case Component::X: U_[cellIdx].setX(eq.phi[cellIdx]); break;
            case Component::Y: U_[cellIdx].setY(eq.phi[cellIdx]); break;
            case Component::Z: U_[cellIdx].setZ(eq.phi[cellIdx]); break;
        }
    }
}
