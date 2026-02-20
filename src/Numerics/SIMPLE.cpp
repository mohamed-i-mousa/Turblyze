
/******************************************************************************
 * @file SIMPLE.cpp
 * @brief Implementation of the SIMPLE algorithm for pressure-velocity coupling
 *****************************************************************************/

#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <map>

#include "LinearInterpolation.hpp"
#include "SIMPLE.hpp"
#include "Scalar.hpp"


// ************************ Constructor & Destructor ************************

SIMPLE::SIMPLE
(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const BoundaryConditions& bc,
    const GradientScheme& gradScheme,
    const ConvectionSchemes& convSchemes
) : allFaces_(faces),
    allCells_(cells),
    bcManager_(bc),
    gradientScheme_(gradScheme),
    convectionScheme_(convSchemes),

    // Physical properties
    rho_(1.225),
    mu_(1.7894e-5),
    nu_(1.7894e-5 / 1.225),

    // Algorithm parameters
    alphaU_(0.7),
    alphaP_(0.3),
    alphaK_(0.5),
    alphaOmega_(0.5),
    maxIterations_(500),
    tolerance_(1e-3),
    enableTurbulence_(false),

    // Turbulence model
    turbulenceModel_(nullptr),
    
    // Field constraint system
    constraintSystem_(nullptr),

    // Solution fields
    U_("U", cells.size(), Vector(0.0, 0.0, 0.0)),
    p_("p", cells.size(), 0.0),
    pCorr_("pCorr", cells.size(), 0.0),
    lastPressureCorrectionRMS_(S(1e9)),

    // Previous-iteration fields
    UPrev_("UPrev", cells.size(), Vector(0.0, 0.0, 0.0)),
    UAvgf_("UAvgf", faces.size(), Vector(0.0, 0.0, 0.0)),
    UAvgPrevf_("UAvgPrevf", faces.size(), Vector(0.0, 0.0, 0.0)),

    // Face-based fields for Rhie-Chow interpolation
    RhieChowFlowRate_("RhieChowMassFlux", faces.size(), 0.0),
    RhieChowFlowRatePrev_("RhieChowMassFluxPrev", faces.size(), 0.0),

    // Momentum equation coefficients
    DU_("DU", cells.size(), 0.0),
    DUf_("DUf", faces.size(), 0.0),

    // Gradient fields
    gradP_("gradP", cells.size(), Vector(0.0, 0.0, 0.0)),
    gradPCorr_("gradPCorr", cells.size(), Vector(0.0, 0.0, 0.0)),

    // Matrix constructor
    matrixConstruct_(nullptr),

    // Linear solvers (per-equation defaults)
    momentumSolver_("momentum", S(1e-8), 1000),
    pressureSolver_("pCorr", S(1e-6), 1000, S(0.05))
{}


// ****************************** Private Helpers ******************************

ScalarField SIMPLE::extractComponent
(
    const VectorField& V,
    int component,
    const std::string& name
)
{
    ScalarField result(name, V.size());

    size_t numCells = V.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        switch (component)
        {
            case 0: result[cellIdx] = V[cellIdx].x(); break;
            case 1: result[cellIdx] = V[cellIdx].y(); break;
            case 2: result[cellIdx] = V[cellIdx].z(); break;
        }
    }

    return result;
}

// ******************************* Public Methods ******************************

void SIMPLE::initialize
(
    const Vector& initialVelocity,
    Scalar initialPressure,
    Scalar initialK,
    Scalar initialOmega
)
{
    matrixConstruct_ = 
        std::make_unique<Matrix>
        (
            allFaces_,
            allCells_,
            bcManager_
        );
        
    // Initialize constraint system
    constraintSystem_ =
        std::make_unique<Constraint>
        (
            U_,
            p_
        );

    U_.setAll(initialVelocity);
    p_.setAll(initialPressure);
    UAvgf_.setAll(initialVelocity);

    // Initialize RhieChowFlowRate_ with linear interpolation
    size_t numFaces = allFaces_.size();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = allFaces_[faceIdx];
        Vector Uf = interpolateToFace(face, U_, bcManager_, "U");
        Vector Sf = face.normal() * face.projectedArea();
        RhieChowFlowRate_[faceIdx] = dot(Uf, Sf);
    }

    if (enableTurbulence_)
    {
        turbulenceModel_ = 
            std::make_unique<kOmegaSST>
            (
                allFaces_,
                allCells_,
                bcManager_,
                gradientScheme_,
                convectionScheme_.k(),
                convectionScheme_.omega()
            );

        turbulenceModel_->
            initialize
            (
                nu_,
                initialK,
                initialOmega,
                alphaK_,
                alphaOmega_
            );

        std::cout
            << "k-omega SST turbulence model initialized." << std::endl;
    }
    
    std::cout
        << "SIMPLE algorithm is initialized!" << std::endl;
}

void SIMPLE::solve()
{
    std::cout
        << "\n=== Starting SIMPLE Loop ===" << std::endl;
    
    int iteration = 0;

    bool converged = false;
    
    while (!converged && iteration < maxIterations_)
    {
        std::cout
            << "\n Iteration " << iteration + 1;

        UPrev_ = U_;
        UAvgPrevf_ = UAvgf_;
        RhieChowFlowRatePrev_ = RhieChowFlowRate_;

        size_t numCells = allCells_.size();
        for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
            gradP_[cellIdx] =
                gradientScheme_.cellGradient(cellIdx, "p", p_);
        }

        solveMomentumEquations();

        calculateRhieChowFlowRate();

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

void SIMPLE::solveMomentumEquations()
{
    ScalarField nuEff("nuEff", allCells_.size());

    // Extract velocity components
    ScalarField Ux = extractComponent(U_, 0, "Ux");
    ScalarField Uy = extractComponent(U_, 1, "Uy");
    ScalarField Uz = extractComponent(U_, 2, "Uz");

    // Extract previous-iteration components
    ScalarField UxPrev = extractComponent(UPrev_, 0, "UxPrev");
    ScalarField UyPrev = extractComponent(UPrev_, 1, "UyPrev");
    ScalarField UzPrev = extractComponent(UPrev_, 2, "UzPrev");

    ScalarField UxSource("UxSource", allCells_.size());
    ScalarField UySource("UySource", allCells_.size());
    ScalarField UzSource("UzSource", allCells_.size());

    // Reset diagonals accumulator
    DU_.setAll(0.0);

    // Build effective viscosity per cell
    size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        if (enableTurbulence_ && turbulenceModel_)
        {
            const ScalarField& nut = turbulenceModel_->getTurbulentViscosity();
            nuEff[cellIdx] = nu_ + nut[cellIdx];
        }
        else
        {
            nuEff[cellIdx] = nu_;
        }

        // Calculate pressure gradients source term
        UxSource[cellIdx] = -gradP_[cellIdx].x() * allCells_[cellIdx].volume();
        UySource[cellIdx] = -gradP_[cellIdx].y() * allCells_[cellIdx].volume();
        UzSource[cellIdx] = -gradP_[cellIdx].z() * allCells_[cellIdx].volume();
    }

    // Reset interpolated diagonals accumulator
    DUf_.setAll(0.0);

    VectorField gradUx("gradUx", numCells, Vector(0.0, 0.0, 0.0));
    VectorField gradUy("gradUy", numCells, Vector(0.0, 0.0, 0.0));
    VectorField gradUz("gradUz", numCells, Vector(0.0, 0.0, 0.0));

    for (size_t cellIdx = 0; cellIdx < allCells_.size(); ++cellIdx)
    {
        gradUx[cellIdx] =
            gradientScheme_.cellGradient(cellIdx, "Ux", Ux);

        gradUy[cellIdx] =
            gradientScheme_.cellGradient(cellIdx, "Uy", Uy);

        gradUz[cellIdx] =
            gradientScheme_.cellGradient(cellIdx, "Uz", Uz);
    }

    // Store velocity gradients for turbulence model
    gradU_ = {gradUx, gradUy, gradUz};

    // Add transpose gradient source term
    if (enableTurbulence_ && turbulenceModel_)
    {
        ScalarField transposeSourceX("transposeSourceX", numCells, 0.0);
        ScalarField transposeSourceY("transposeSourceY", numCells, 0.0);
        ScalarField transposeSourceZ("transposeSourceZ", numCells, 0.0);

        calculateTransposeGradientSource
        (
            nuEff,
            gradUx,
            gradUy,
            gradUz,
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
    solveMomentumComponent('x', Ux, UxSource, UxPrev, nuEff, gradUx);
    solveMomentumComponent('y', Uy, UySource, UyPrev, nuEff, gradUy);
    solveMomentumComponent('z', Uz, UzSource, UzPrev, nuEff, gradUz);

    // Calculate DUf_ field using complete DU_
    size_t numFaces = allFaces_.size();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = allFaces_[faceIdx];

        if (face.isBoundary())
        {
            const BoundaryData* bc =
                bcManager_.fieldBC(face.patch()->patchName(), "p");

            if (bc && bc->type() == BCType::FIXED_VALUE)
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

void SIMPLE::solveMomentumComponent
(
    char component,
    ScalarField& UComponent,
    const ScalarField& source,
    const ScalarField& UComponentPrev,
    const ScalarField& nuEff,
    const VectorField& gradComponent
)
{
    std::string fieldName = "U";
    fieldName += component;

    matrixConstruct_->
        buildMatrix
        (
            UComponent,
            source,
            RhieChowFlowRatePrev_,
            nuEff,
            convectionScheme_.momentum(),
            gradientScheme_,
            gradComponent,
            fieldName
        );

    // Get references to the assembled matrix
    auto& matrixA = matrixConstruct_->getMatrixA();
    auto& vectorB = matrixConstruct_->getVectorB();

    matrixConstruct_->relax(alphaU_, UComponentPrev);

    // Store DU_
    size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        DU_[cellIdx] = 
            allCells_[cellIdx].volume()
          / (matrixA.coeff(cellIdx, cellIdx) + vSmallValue);
    }

    // Map UComponent field storage as Eigen vector (zero-copy)
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    solution(UComponent.data(), numCells);

    momentumSolver_.solveWithBiCGSTAB
    (
        solution,
        matrixA,
        vectorB
    );

    // Write solved component back to velocity VectorField
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        switch (component)
        {
            case 'x': U_[cellIdx].setX(UComponent[cellIdx]); break;
            case 'y': U_[cellIdx].setY(UComponent[cellIdx]); break;
            case 'z': U_[cellIdx].setZ(UComponent[cellIdx]); break;
        }
    }
}

void SIMPLE::calculateRhieChowFlowRate()
{
    size_t numFaces = allFaces_.size();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = allFaces_[faceIdx];

        if (face.isBoundary())
        {
            UAvgf_[faceIdx] =
                bcManager_.calculateBoundaryVectorFaceValue(face, U_, "U");

            RhieChowFlowRate_[faceIdx] =
                dot(UAvgf_[faceIdx], face.normal() * face.projectedArea());

            continue;
        }

        const size_t P = face.ownerCell();
        const size_t N = face.neighborCell().value();

        // Linear-interpolated velocity at face
        const Vector UfLinear =
            interpolateToFace(face, U_, bcManager_, "U");

        const Vector gradPAvgf = 
            interpolateToFace(face, gradP_);
                
        const Vector Sf = face.normal() * face.projectedArea();
        
        Vector gradPf =
            gradientScheme_.faceGradient
            (
                faceIdx,
                "p",
                p_,
                gradP_[P],
                gradP_[N]
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
    size_t numCells = allCells_.size();
    VectorField gradPCorrPrecomputed
    (
        "gradPCorr", numCells, Vector(0.0, 0.0, 0.0)
    );

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradPCorrPrecomputed[cellIdx] =
            gradientScheme_.cellGradient(cellIdx, "pCorr", pCorr_);
    }

    matrixConstruct_->
        buildPressureCorrectionMatrix
        (
            RhieChowFlowRate_,
            DUf_,
            pCorr_,
            gradientScheme_,
            gradPCorrPrecomputed
        );

    // Get references to the assembled matrix
    auto& matrixA = matrixConstruct_->getMatrixA();
    auto& vectorB = matrixConstruct_->getVectorB();

    // Map pCorr field storage as Eigen vector (zero-copy)
    pCorr_.setAll(S(0.0));
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>
    pCorrSolution(pCorr_.data(), numCells);

    pressureSolver_.solveWithPCG
    (
        pCorrSolution,
        matrixA,
        vectorB
    );

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradPCorr_[cellIdx] =
            gradientScheme_.cellGradient(cellIdx, "pCorr", pCorr_);
    }
}

void SIMPLE::correctVelocity()
{
    size_t numCells = allCells_.size();
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
    
    if (constraintSystem_)
    {
        int velocityConstraints = 
            constraintSystem_->applyVelocityConstraints();

        if (velocityConstraints > 0)
        {
            std::cout
                << "  Applied velocity constraints to "
                << velocityConstraints << " cells" << std::endl;
        }
    }

    // Update face velocities
    size_t numFaces = allFaces_.size();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = allFaces_[faceIdx];

        if (face.isBoundary())
        {
            UAvgf_[faceIdx] =
                bcManager_.calculateBoundaryVectorFaceValue(face, U_, "U");
        }
        else
        {
            UAvgf_[faceIdx] =
                interpolateToFace(face, U_, bcManager_, "U");
        }
    }
}

void SIMPLE::correctPressure()
{
    Scalar sumSq = 0.0;

    size_t numCells = allCells_.size();
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx) 
    {
        sumSq += pCorr_[cellIdx] * pCorr_[cellIdx];
    }

    lastPressureCorrectionRMS_ = std::sqrt(sumSq / numCells);

    // Apply pressure correction
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        p_[cellIdx] += alphaP_ * pCorr_[cellIdx];
    }
    
    // Apply pressure bounds constraints
    if (constraintSystem_)
    {
        int pressureConstraints = 
            constraintSystem_->applyPressureConstraints();

        if (pressureConstraints > 0)
        {
            std::cout
                << "  Applied pressure constraints to "
                << pressureConstraints << " cells" << std::endl;
        }
    }

    // Reset pressure correction for next iteration
    pCorr_.setAll(0.0);
}

void SIMPLE::correctFlowRate()
{
    Scalar flowRateCorrection;

    // Update mass flux on faces
    size_t numFaces = allFaces_.size();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = allFaces_[faceIdx];

        if (face.isBoundary())
        {
            const BoundaryData* bc =
                bcManager_.fieldBC(face.patch()->patchName(), "p");

            if (bc && (bc->type() == BCType::FIXED_VALUE
                    || bc->type() == BCType::ZERO_GRADIENT))
            {
                continue;
            }

            Scalar gradn =
                dot(gradPCorr_[face.ownerCell()], face.normal());

            flowRateCorrection =
                DU_[face.ownerCell()] * gradn * face.projectedArea();

            RhieChowFlowRate_[faceIdx] -= flowRateCorrection;
            continue;
        }

        size_t ownerIdx = face.ownerCell();
        size_t neighborIdx = face.neighborCell().value();

        Vector gradPCorrf =
            gradientScheme_.faceGradient
            (
                faceIdx,
                "pCorr",
                pCorr_,
                gradPCorr_[ownerIdx],
                gradPCorr_[neighborIdx]
            );

        Vector Sf = face.normal() * face.projectedArea();

        flowRateCorrection = DUf_[faceIdx] * dot(gradPCorrf, Sf);

        RhieChowFlowRate_[faceIdx] -= flowRateCorrection;
    }
}

void SIMPLE::solveTurbulence()
{
    if (enableTurbulence_ && turbulenceModel_)
    {
        std::cout
            << "  Solving turbulence equations..." << std::endl;

        // Snapshot k/omega before solve for residual tracking
        const auto& kField = turbulenceModel_->getk();
        const auto& omegaField =
            turbulenceModel_->getOmega();
        size_t numCells = allCells_.size();

        std::vector<Scalar> kPrev(numCells);
        std::vector<Scalar> omegaPrev(numCells);
        for (size_t i = 0; i < numCells; ++i)
        {
            kPrev[i] = kField[i];
            omegaPrev[i] = omegaField[i];
        }

        turbulenceModel_->
            solve
            (
                U_,
                RhieChowFlowRate_,
                gradU_
            );

        // Compute normalised change: ||x - x_prev|| / ||x_prev||
        Scalar kDiffSq = 0.0, kPrevSq = 0.0;
        Scalar omDiffSq = 0.0, omPrevSq = 0.0;
        for (size_t i = 0; i < numCells; ++i)
        {
            Scalar dk = kField[i] - kPrev[i];
            kDiffSq += dk * dk;
            kPrevSq += kPrev[i] * kPrev[i];

            Scalar dw = omegaField[i] - omegaPrev[i];
            omDiffSq += dw * dw;
            omPrevSq += omegaPrev[i] * omegaPrev[i];
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
    // Check mass conservation and velocity/pressure residuals
    Scalar massImbalance = calculateMassImbalance();
    Scalar velocityResidual = calculateVelocityResidual();
    Scalar pressureResidual = calculatePressureResidual();

    bool converged =
        (
            (massImbalance < tolerance_)
         && (velocityResidual < tolerance_)
         && (pressureResidual < tolerance_)
        );

    // Include turbulence residuals when enabled
    if (enableTurbulence_)
    {
        converged = converged
            && (lastKResidual_ < tolerance_)
            && (lastOmegaResidual_ < tolerance_);
    }

    std::cout
        << " - Mass: " << std::scientific << massImbalance
        << ", Velocity: " << velocityResidual
        << ", Pressure: " << pressureResidual;

    if (enableTurbulence_)
    {
        std::cout
            << ", k: " << lastKResidual_
            << ", omega: " << lastOmegaResidual_;
    }

    std::cout
        << std::fixed << std::endl;

    return converged;
}


// ****************************** Private Methods ******************************

Scalar SIMPLE::calculateMassImbalance() const
{
    // Dimensionless normalized continuity residual per cell, averaged
    Scalar totalNormImbalance = 0.0;

    size_t numCells = allCells_.size();
    for (size_t i = 0; i < numCells; ++i) 
    {
        Scalar net = 0.0;
        Scalar sumAbs = 0.0;

        for (size_t j = 0; j < allCells_[i].faceIndices().size(); ++j)
        {
            const size_t faceIdx = allCells_[i].faceIndices()[j];
            const int sign = allCells_[i].faceSigns()[j];
            const Scalar mf = RhieChowFlowRate_[faceIdx];
            net += sign * mf;
            sumAbs += std::abs(mf);
        }

        const Scalar denom = sumAbs + vSmallValue;
        totalNormImbalance += std::abs(net) / denom;
    }
    
    return totalNormImbalance / std::max<size_t>(1, numCells);
}

Scalar SIMPLE::calculateVelocityResidual() const
{
    // Normalized residual: ||U - U_prev||_2 / (||U_prev||_2 + eps)
    Scalar num = 0.0;
    Scalar den = 0.0;

    for (size_t i = 0; i < allCells_.size(); ++i)
    {
        const Vector d = U_[i] - UPrev_[i];
        num += d.magnitudeSquared();
        den += UPrev_[i].magnitudeSquared();
    }

    num = std::sqrt(num + vSmallValue);
    den = std::sqrt(den + vSmallValue);

    return num / den;
}

Scalar SIMPLE::calculatePressureResidual() const
{
    // Normalize p' RMS by RMS(p)
    Scalar sumP2 = 0.0;

    for (size_t i = 0; i < allCells_.size(); ++i) 
    {
        sumP2 += p_[i] * p_[i];
    }

    const Scalar pRms = std::sqrt(sumP2 / (allCells_.size()));

    return lastPressureCorrectionRMS_ / (pRms + vSmallValue);
}

void SIMPLE::calculateTransposeGradientSource
(
    const ScalarField& nuEff,
    const VectorField& gradUx,
    const VectorField& gradUy,
    const VectorField& gradUz,
    ScalarField& transposeSourceX,
    ScalarField& transposeSourceY,
    ScalarField& transposeSourceZ
)
{
    transposeSourceX.setAll(0.0);
    transposeSourceY.setAll(0.0);
    transposeSourceZ.setAll(0.0);

    for (size_t faceIdx = 0; faceIdx < allFaces_.size(); ++faceIdx)
    {
        const Face& face = allFaces_[faceIdx];
        const size_t ownerIdx = face.ownerCell();
        Vector Sf = face.normal() * face.projectedArea();

        if (face.isBoundary())
        {
            // For boundary faces, use owner cell gradients
            Scalar nuEfff = nuEff[ownerIdx];

            // Build transpose gradient columns
            Vector gradTx = 
                buildGradientTransposeColumn
                (
                    gradUx[ownerIdx], 
                    gradUy[ownerIdx], 
                    gradUz[ownerIdx], 
                    0
                );

            Vector gradTy = 
                buildGradientTransposeColumn
                (
                    gradUx[ownerIdx], 
                    gradUy[ownerIdx], 
                    gradUz[ownerIdx], 
                    1
                );

            Vector gradTz = 
                buildGradientTransposeColumn
                (
                    gradUx[ownerIdx], 
                    gradUy[ownerIdx], 
                    gradUz[ownerIdx], 
                    2
                );

            // Flux contribution: ν_eff_f * (∇U)^T · Sf
            transposeSourceX[ownerIdx] += nuEfff * dot(gradTx, Sf);
            transposeSourceY[ownerIdx] += nuEfff * dot(gradTy, Sf);
            transposeSourceZ[ownerIdx] += nuEfff * dot(gradTz, Sf);
        }
        else
        {
            // Internal face: interpolate gradients to face
            const size_t neighborIdx = face.neighborCell().value();

            // Interpolate effective viscosity and gradients to face
            Scalar nuEfff = interpolateToFace(face, nuEff);
            Vector gradUxf = interpolateToFace(face, gradUx);
            Vector gradUyf = interpolateToFace(face, gradUy);
            Vector gradUzf = interpolateToFace(face, gradUz);

            // Build transpose gradient columns
            Vector gradTx =    
                buildGradientTransposeColumn
                (
                    gradUxf, 
                    gradUyf, 
                    gradUzf, 
                    0
                );

            Vector gradTy = 
                buildGradientTransposeColumn
                (
                    gradUxf, 
                    gradUyf, 
                    gradUzf, 
                    1
                );

            Vector gradTz = 
                buildGradientTransposeColumn
                (
                    gradUxf, 
                    gradUyf, 
                    gradUzf, 
                    2
                );

            // Flux contribution: ν_eff_f * (∇U)^T · Sf
            Scalar fluxX = nuEfff * dot(gradTx, Sf);
            Scalar fluxY = nuEfff * dot(gradTy, Sf);
            Scalar fluxZ = nuEfff * dot(gradTz, Sf);

            transposeSourceX[ownerIdx] += fluxX;
            transposeSourceY[ownerIdx] += fluxY;
            transposeSourceZ[ownerIdx] += fluxZ;

            transposeSourceX[neighborIdx] -= fluxX;
            transposeSourceY[neighborIdx] -= fluxY;
            transposeSourceZ[neighborIdx] -= fluxZ;
        }
    }
}

Vector SIMPLE::buildGradientTransposeColumn
(
    const Vector& gradUx,
    const Vector& gradUy,
    const Vector& gradUz,
    int component
) const
{
    switch (component)
    {
        case 0: return Vector(gradUx.x(), gradUy.x(), gradUz.x());
        case 1: return Vector(gradUx.y(), gradUy.y(), gradUz.y());
        case 2: return Vector(gradUx.z(), gradUy.z(), gradUz.z());
        default: return Vector(0.0, 0.0, 0.0);
    }
}