
/******************************************************************************
 * @file SIMPLE.cpp
 * @brief Implementation of the SIMPLE algorithm for pressure-velocity coupling
 *****************************************************************************/

#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <map>

#include "linearInterpolation.h"
#include "SIMPLE.h"
#include "Scalar.h"
#include "flowRate.h"

SIMPLE::SIMPLE
(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells, 
    const BoundaryConditions& bc,
    const GradientScheme& gradScheme,
    const ConvectionScheme& convScheme
) : allFaces(faces),
    allCells(cells),
    bcManager(bc),
    gradientScheme(gradScheme),
    convectionScheme(convScheme),

    // Physical properties
    rho(1.225),
    mu(1.7894e-5),
    nu(1.7894e-5 / 1.225),

    // Algorithm parameters
    alpha_U(0.7),
    alpha_p(0.3),
    maxIterations(500),
    tolerance(1e-3),
    enableTurbulence(false),

    // Turbulence model
    turbulenceModel(nullptr),
    
    // Field constraint system
    constraintSystem(nullptr),

    // Solution fields
    U("U", cells.size(), Vector(0.0, 0.0, 0.0)),
    p("p", cells.size(), 0.0),
    pCorr("pCorr", cells.size(), 0.0),
    lastPressureCorrectionRMS(S(1e9)),

    // Previous-iteration fields
    U_prev("U_prev", cells.size(), Vector(0.0, 0.0, 0.0)),
    U_x_prev("U_x_prev", cells.size(), 0.0),
    U_y_prev("U_y_prev", cells.size(), 0.0),
    U_z_prev("U_z_prev", cells.size(), 0.0),
    U_f_avg("U_face", faces.size(), Vector(0.0, 0.0, 0.0)),
    U_f_avg_prev("U_face_prev", faces.size(), Vector(0.0, 0.0, 0.0)),

    // Face-based fields for Rhie-Chow interpolation
    RhieChowFlowRate("RhieChowMassFlux", faces.size(), 0.0),
    RhieChowFlowRate_prev("RhieChowMassFlux_prev", faces.size(), 0.0),

    // Momentum equation coefficients
    DU("DU", cells.size(), 0.0),
    DUf("DUf", faces.size(), 0.0),

    // Gradient fields
    gradP("gradP", cells.size(), Vector(0.0, 0.0, 0.0)),
    gradPCorr("gradPCorr", cells.size(), Vector(0.0, 0.0, 0.0)),

    // Matrix constructor
    matrixConstruct(nullptr)
{
    initialize(Vector(0.0, 0.0, -0.5), 0.0);
}

void SIMPLE::initialize(const Vector& initialVelocity, Scalar initialPressure)
{
    matrixConstruct = 
        std::make_unique<Matrix>
        (
            allFaces,
            allCells,
            bcManager
        );
        
    // Initialize constraint system
    constraintSystem = std::make_unique<Constraint>
    (
        U,        // velocity field
        p         // pressure field
    );
    
    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        U[cellIdx] = initialVelocity;
        p[cellIdx] = initialPressure;                     
    }

    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
    {
        U_f_avg[faceIdx] = initialVelocity;
    }

    if (enableTurbulence)
    {
        turbulenceModel = std::make_unique<KOmegaSST>
        (
            allFaces,
            allCells,
            bcManager,
            gradientScheme
        );

        turbulenceModel->initialize(U, nu);

        std::cout << "k-omega SST turbulence model initialized." << std::endl;
    }
    
    std::cout << "SIMPLE algorithm is initialized!" << std::endl;
}

void SIMPLE::solve()
{
    std::cout << "\n=== Starting SIMPLE Loop ===" << std::endl;
    
    int iteration = 0;

    bool converged = false;
    
    while (!converged && iteration < maxIterations)
    {
        std::cout   << "\n--- SIMPLE Iteration " << iteration + 1
                    << " ---" << std::endl;

        U_prev = U;
        U_f_avg_prev = U_f_avg;
        RhieChowFlowRate_prev = RhieChowFlowRate;

        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            // Compute pressure gradients for this iteration
            gradP[cellIdx] = gradientScheme.CellGradient(cellIdx, p, allCells);

            // Get previous iteration value for velocity components
            U_x_prev[cellIdx] = U_prev[cellIdx].x();
            U_y_prev[cellIdx] = U_prev[cellIdx].y();
            U_z_prev[cellIdx] = U_prev[cellIdx].z();
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
        std::cout   << "WARNING: SIMPLE algorithm did not converge after "
                    << maxIterations << " iterations." << std::endl;
    }
    else
    {
        std::cout   << "SIMPLE algorithm converged in " << iteration
                    << " iterations." << std::endl;
    }
}

void SIMPLE::solveMomentumEquations() 
{
    ScalarField nu_eff("nu_eff", allCells.size());

    ScalarField U_x("U_x", allCells.size());
    ScalarField U_y("U_y", allCells.size());
    ScalarField U_z("U_z", allCells.size());

    ScalarField Ux_source("U_x_source", allCells.size());
    ScalarField Uy_source("U_y_source", allCells.size());
    ScalarField Uz_source("U_z_source", allCells.size());

    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        // Reset diagonals accumulator 
        DU[cellIdx] = 0.0;

        // Build effictive viscosity per cell. μ_eff = nu + nu_t
        if (enableTurbulence && turbulenceModel) 
        {
            const ScalarField& nu_t = turbulenceModel->getTurbulentViscosity();
            nu_eff[cellIdx] = nu + nu_t[cellIdx];
        }
        else
        {
            nu_eff[cellIdx] = nu;
        }

        // Extract momentum components
        U_x[cellIdx] = U[cellIdx].x();
        U_y[cellIdx] = U[cellIdx].y();
        U_z[cellIdx] = U[cellIdx].z();

        // Calculate pressure gradients source term: S_i = -(∂p/∂x_i) * V
        Ux_source[cellIdx] = -gradP[cellIdx].x() * allCells[cellIdx].volume();
        Uy_source[cellIdx] = -gradP[cellIdx].y() * allCells[cellIdx].volume();
        Uz_source[cellIdx] = -gradP[cellIdx].z() * allCells[cellIdx].volume();
    }
    
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
    {
        // Reset interpolated diagonals accumulator
        DUf[faceIdx] = 0.0;
    }

    // ****************************** x-momentum ******************************

    matrixConstruct->buildMatrix
    (
        U_x,
        Ux_source,
        U_prev,
        nu_eff,
        convectionScheme,
        gradientScheme,
        "U_x"
    );
    
    // Get references to the assembled matrix and vector
    auto& b_vector = 
        const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>
        (
            matrixConstruct->getVectorB()
        );

    auto& A_matrix = 
        const_cast<Eigen::SparseMatrix<Scalar>&>
        (
            matrixConstruct->getMatrixA()
        );

    // Store the diagonal coefficient for Rhie-Chow interpolation
    for (size_t i = 0; i < allCells.size(); ++i)
    {
        DU[i] = allCells[i].volume() / (A_matrix.coeff(i, i) + vSmallValue);
    }

    // Apply implicit under-relaxation before solving (Patankar's relaxation)
    matrixConstruct->relax(alpha_U, U_x_prev);

    // Solve for U.x
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> U_x_solution(allCells.size());

    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        U_x_solution(cellIdx) = U[cellIdx].x();
    }
    
    bool solved = LinearSolvers::BiCGSTAB
    (
        U_x_solution, 
        A_matrix, 
        b_vector, 
        1e-8, 
        1000, 
        "U_x"
    );
    
    if (solved)
    {
        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx) 
        {
            U[cellIdx].setX(U_x_solution(cellIdx));
        }
    }
    
    // ****************************** y-momentum ******************************

    // Build y-momentum matrix
    matrixConstruct->buildMatrix
    (
        U_y,
        Uy_source,
        U_prev,
        nu_eff,
        convectionScheme,
        gradientScheme,
        "U_y"
    );
    
    // Get references to the assembled matrix and vector
    auto& b_vector_y = 
        const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>
        (
            matrixConstruct->getVectorB()
        );

    auto& A_matrix_y = 
        const_cast<Eigen::SparseMatrix<Scalar>&>
        (
            matrixConstruct->getMatrixA()
        );

    // Apply implicit under-relaxation before solving (Patankar's relaxation)
    matrixConstruct->relax(alpha_U, U_y_prev);
    
    // Solve for U.y
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> U_y_solution(allCells.size());
    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        U_y_solution(cellIdx) = U[cellIdx].y();
    }
    
    solved = LinearSolvers::BiCGSTAB
    (
        U_y_solution,
        A_matrix_y,
        b_vector_y,
        1e-8,
        1000,
        "U_y"
    );
    
    if (solved)
    {
        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            U[cellIdx].setY(U_y_solution(cellIdx));
        }
    }
    
    // ****************************** z-momentum ******************************

    // Build z-momentum matrix
    matrixConstruct->buildMatrix
    (
        U_z,
        Uz_source,
        U_prev,
        nu_eff,
        convectionScheme,
        gradientScheme,
        "U_z"
    );
    
    // Get references to the assembled matrix and vector
    auto& b_vector_z = 
        const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>
        (
            matrixConstruct->getVectorB()
        );

    auto& A_matrix_z = 
        const_cast<Eigen::SparseMatrix<Scalar>&>
        (
           matrixConstruct->getMatrixA()
        );

    // Apply implicit under-relaxation before solving (Patankar's relaxation)
    matrixConstruct->relax(alpha_U, U_z_prev);

    // Solve for U.z
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> U_z_solution(allCells.size());

    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        U_z_solution(cellIdx) = U[cellIdx].z();
    }
    
    solved = LinearSolvers::BiCGSTAB
    (
        U_z_solution,
        A_matrix_z,
        b_vector_z,
        1e-8,
        1000,
        "U_z"
    );
    
    if (solved)
    {
        for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
        {
            U[cellIdx].setZ(U_z_solution(cellIdx));
        }
    }
    
    // Calculate DUf field using complete DU
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
    {
        const Face& face = allFaces[faceIdx];
        
        if (face.isBoundary())
        {
            // Find which patch this face belongs to
            std::string patchName = "";
            for (const auto& patch : bcManager.patches())
            {
                if 
                (
                    face.id() >= patch.firstFaceIndex() 
                 && face.id() <= patch.lastFaceIndex()
                )
                {
                    patchName = patch.patchName();
                    break;
                }
            }
            
            const BoundaryData* bc = bcManager.fieldBC(patchName, "p");
            
            if (bc && bc->type() == BCType::FIXED_VALUE)
            {
                // Fixed pressure boundary: normal pressure-velocity coupling
                DUf[faceIdx] = DU[allFaces[faceIdx].ownerCell()];
            }
            else
            {
                // Zero gradient pressure boundary: no pressure-velocity coupling
                DUf[faceIdx] = 0.0;
            }
        }
        else
        {
            // Linear interpolation of momentum diagonal coefficients to face
            DUf[faceIdx] = linearInterpolation(face, DU);
        }
    }
}

void SIMPLE::calculateRhieChowFlowRate()
{
    // U_f = U_f_linearInterpolated
    //     + DUf * (∇p_face_linearInterpolated - ∇p_f)
    //     + (1 - alpha_U) * (U_f_previousIter - U_f_previousIter_linearInterpolated)

    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
    {
        const Face& face = allFaces[faceIdx];

        if (face.isBoundary()) 
        {
            // Centralized boundary condition handling for velocity
            U_f_avg[faceIdx] = 
                bcManager.calculateBoundaryFaceVectorValue(face, U, "U");

            RhieChowFlowRate[faceIdx] =
                dot(U_f_avg[faceIdx], face.normal() * face.area());
                
            continue;
        }

        const size_t P = face.ownerCell();
        const size_t N = face.neighborCell().value();

        // Linear-interpolated values with proper boundary condition handling
        const Vector U_f_linear = 
            VectorLinearInterpolation(face, U, bcManager, "U");

        const Vector gradP_f_avg = 
            VectorLinearInterpolation(face, gradP);
                
        const Vector S_f = face.normal() * face.area();
        
        Vector gradP_f = 
            gradientScheme.FaceGradient
            (
                faceIdx,
                gradP[P],
                gradP[N],
                p,
                allCells,
                allFaces,
                bcManager,
                "p"
            );

        RhieChowFlowRate[faceIdx] = 
            dot(U_f_linear, S_f) 
          - dot((DUf[faceIdx] * (gradP_f - gradP_f_avg)), S_f)
          + (S(1.0) - alpha_U) 
          * (RhieChowFlowRate_prev[faceIdx] - dot(U_f_avg_prev[faceIdx], S_f));
    }
}

void SIMPLE::solvePressureCorrection()
{
    // Build pressure correction equation
    // ∇·(DUf ∇p') = ∇·ρU  >> ∇·(DUf ∇p') = -∑ m*_f
    matrixConstruct->buildPressureCorrectionMatrix
    (
        RhieChowFlowRate,
        DUf
    );

    auto& A_matrix = const_cast<Eigen::SparseMatrix<Scalar>&>
    (
        matrixConstruct->getMatrixA()
    );

    auto& b_vector = const_cast<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&>
    (
        matrixConstruct->getVectorB()
    );

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> pCorr_solution(allCells.size());

    pCorr_solution.setZero();

    bool solved = LinearSolvers::PCG
    (
        pCorr_solution,
        A_matrix,
        b_vector,
        1e-8,
        1000,
        "pCorr"
    );

    if (solved)
    {
        for (size_t i = 0; i < allCells.size(); ++i)
        {
            pCorr[i] = pCorr_solution(i);
        }
        
        // Compute pressure correction gradients for velocity correction
        for (size_t i = 0; i < allCells.size(); ++i)
        {
            gradPCorr[i] = gradientScheme.CellGradient(i, pCorr, allCells);
        }
    }
}

void SIMPLE::correctVelocity()
{
    // Velocity correction: U = U* - DU * ∇p'
    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        U[cellIdx].setX(U[cellIdx].x() - DU[cellIdx] * gradPCorr[cellIdx].x());
        U[cellIdx].setY(U[cellIdx].y() - DU[cellIdx] * gradPCorr[cellIdx].y());
        U[cellIdx].setZ(U[cellIdx].z() - DU[cellIdx] * gradPCorr[cellIdx].z());
    }
    
    // Apply velocity field constraints after correction
    if (constraintSystem)
    {
        int velocityConstraints = constraintSystem->applyVelocityConstraints();
        if (velocityConstraints > 0)
        {
            std::cout   << "  Applied velocity constraints to " 
                        << velocityConstraints << " cells" << std::endl;
        }
    }
    
    // Update face velocities to match corrected cell velocities
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
    {
        const Face& face = allFaces[faceIdx];
        if (face.isBoundary()) 
        {
            U_f_avg[faceIdx] = 
                bcManager.calculateBoundaryFaceVectorValue(face, U, "U");
        }
        else
        {
            U_f_avg[faceIdx] = 
                VectorLinearInterpolation(face, U, bcManager, "U");
        }
    }
}

void SIMPLE::correctPressure()
{
    // Pressure correction: p = p + α_p * p'
    // Track p' RMS before reset to report meaningful pressure residual
    Scalar sumSq = 0.0;

    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx) 
    {
        sumSq += pCorr[cellIdx] * pCorr[cellIdx];
    }

    lastPressureCorrectionRMS = std::sqrt(sumSq / allCells.size());

    // Apply pressure correction
    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        p[cellIdx] += alpha_p * pCorr[cellIdx];
    }
    
    // Apply pressure bounds constraints if enabled
    if (constraintSystem)
    {
        int pressureConstraints = 
            constraintSystem->applyPressureConstraints();

        if (pressureConstraints > 0)
        {
            std::cout   << "  Applied pressure constraints to " 
                        << pressureConstraints << " cells" << std::endl;
        }
    }

    // Reset pressure correction for next iteration
    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx)
    {
        pCorr[cellIdx] = 0.0;
    }
}

void SIMPLE::correctFlowRate() 
{       
    Scalar flowRateCorrection; 

    // Update mass flux on faces to be consistent with corrected pressure
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx) 
    {
        const Face& face = allFaces[faceIdx];

        if (face.isBoundary()) 
        {
            std::string patchName = "";
            for (const auto& patch : bcManager.patches()) 
            {
                if (face.id() >= patch.firstFaceIndex() 
                 && face.id() <= patch.lastFaceIndex()) 
                    {
                    patchName = patch.patchName();
                    break;
                }
            }
        
            const BoundaryData* bc = bcManager.fieldBC(patchName, "p");
            
            if (bc && bc->type() == BCType::FIXED_VALUE)
            {
                // Skip fixed pressure boundaries
                continue;
            }
            
            // Apply pressure correction for non-fixed pressure boundaries
            Vector grad_pCorr = 
                gradientScheme.CellGradient(face.ownerCell(), pCorr, allCells);
            
            Scalar normal_grad = dot(grad_pCorr, face.normal());

            flowRateCorrection = 
                alpha_p * DU[face.ownerCell()] * normal_grad * face.area();
            
            RhieChowFlowRate[faceIdx] -= flowRateCorrection;

            continue;
        }

        // Calculate face gradient of pressure correction
        size_t ownerIdx = face.ownerCell();
        size_t neighborIdx = face.neighborCell().value();
        
        // Face gradient: (pCorr_N - pCorr_P) / |d_PN|
        Vector d_PN = 
            allCells[neighborIdx].centroid() - allCells[ownerIdx].centroid();

        Scalar d_PN_mag = d_PN.magnitude();

        Scalar gradpCorr_f = 
            (pCorr[neighborIdx] - pCorr[ownerIdx]) / (d_PN_mag + vSmallValue);
            
        flowRateCorrection = 
            alpha_p * DUf[faceIdx] * gradpCorr_f * face.area();

        RhieChowFlowRate[faceIdx] -= flowRateCorrection;
    }
}

void SIMPLE::solveTurbulence()
{
    if (enableTurbulence && turbulenceModel)
    {
        std::cout << "  Solving turbulence equations..." << std::endl;
        
        turbulenceModel->solve(U, nu);
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
            (massImbalance < tolerance) 
         && (velocityResidual < tolerance) 
         && (pressureResidual < tolerance)
        );
    
    std::cout   << "  Residuals - Mass: " << std::scientific << massImbalance
                << ", Velocity: " << velocityResidual
                << ", Pressure: " << pressureResidual
                << std::fixed << std::endl;
    
    // Additional convergence monitoring
    if 
    (
        massImbalance > 1e3 
     || velocityResidual > 1e3 
     || pressureResidual > 1e3
    )
    {
        std::cout   << "  WARNING: Residuals are very large -"
                    << " solution may be diverging!" << std::endl;

        std::cout   << "  Consider: " << std::endl;

        std::cout   << "    - Reducing relaxation factors (current: U="
                    << alpha_U << ", p=" << alpha_p << ")" << std::endl;

        std::cout   << "    - Checking boundary conditions" << std::endl;

        std::cout   << "    - Improving initial guess" << std::endl;
    }
    
    return converged;
}

Scalar SIMPLE::calculateMassImbalance() const
{
    // Dimensionless normalized continuity residual per cell, averaged
    Scalar totalNormImbalance = 0.0;

    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        Scalar net = 0.0;
        Scalar sumAbs = 0.0;

        for (size_t j = 0; j < allCells[i].faceIndices().size(); ++j)
        {
            const size_t faceIdx = allCells[i].faceIndices()[j];
            const int sign = allCells[i].faceSigns()[j];
            const Scalar mf = RhieChowFlowRate[faceIdx];
            net += sign * mf;
            sumAbs += std::abs(mf);
        }

        const Scalar denom = sumAbs + vSmallValue;
        totalNormImbalance += std::abs(net) / denom;
    }
    
    return totalNormImbalance / std::max<size_t>(1, allCells.size());
}

Scalar SIMPLE::calculateVelocityResidual() const
{
    // Normalized delta-based residual: ||U - U_prev||_2 / (||U_prev||_2 + eps)
    Scalar num = 0.0;
    Scalar den = 0.0;

    for (size_t i = 0; i < allCells.size(); ++i)
    {
        const Vector d = U[i] - U_prev[i];
        num += d.magnitudeSquared();
        den += U_prev[i].magnitudeSquared();
    }

    num = std::sqrt(num + vSmallValue);
    den = std::sqrt(den + vSmallValue);

    return num / den;
}

Scalar SIMPLE::calculatePressureResidual() const
{
    // Normalize p' RMS by RMS(p)
    Scalar sumP2 = 0.0;

    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        sumP2 += p[i] * p[i];
    }

    const Scalar p_rms = std::sqrt(sumP2 / (allCells.size()));

    return lastPressureCorrectionRMS / (p_rms + vSmallValue);
}

// Setter methods
void SIMPLE::setRelaxationFactors(Scalar alpha_U_new, Scalar alpha_p_new)
{
    alpha_U = alpha_U_new;
    alpha_p = alpha_p_new;
}

void SIMPLE::setConvergenceTolerance(Scalar tol) 
{
    tolerance = tol;
}

void SIMPLE::setMaxIterations(int maxIter) 
{
    maxIterations = maxIter;
}

void SIMPLE::enableTurbulenceModeling(bool enable) 
{
    enableTurbulence = enable;

    if (enable) 
    {
        std::cout   << "k-omega SST turbulence modeling enabled." 
                    << std::endl;
    } 
    else 
    {
        std::cout   << "Laminar flow (turbulence modeling disabled)." 
                    << std::endl;
    }
}

// Turbulence getters
const ScalarField* SIMPLE::getTurbulentKineticEnergy() const 
{
    if (enableTurbulence && turbulenceModel)
    {
        return &(turbulenceModel->getK());
    }

    return nullptr;
}

const ScalarField* SIMPLE::getSpecificDissipationRate() const
{
    if (enableTurbulence && turbulenceModel)
    {
        return &(turbulenceModel->getOmega());
    }

    return nullptr;
}

const ScalarField* SIMPLE::getTurbulentViscosity() const
{
    if (enableTurbulence && turbulenceModel)
    {
        return &(turbulenceModel->getTurbulentViscosity());
    }

    return nullptr;
}

const ScalarField* SIMPLE::getWallDistance() const 
{
    if (enableTurbulence && turbulenceModel)
    {
        return &(turbulenceModel->getWallDistance());
    }

    return nullptr;
}

// Velocity component getters
ScalarField SIMPLE::getVelocityX() const 
{
    ScalarField U_x("U_x", allCells.size());

    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        U_x[i] = U[i].x();
    }

    return U_x;
}

ScalarField SIMPLE::getVelocityY() const 
{
    ScalarField U_y("U_y", allCells.size());

    for (size_t i = 0; i < allCells.size(); ++i) 
    {
        U_y[i] = U[i].y();
    }

    return U_y;
}

ScalarField SIMPLE::getVelocityZ() const 
{
    ScalarField U_z("U_z", allCells.size());

    for (size_t i = 0; i < allCells.size(); ++i)
    {
        U_z[i] = U[i].z();
    }
    
    return U_z;
}


