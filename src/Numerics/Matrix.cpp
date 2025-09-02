#include "Matrix.h"
#include "KOmegaSST.h"
#include <iostream>
#include <algorithm>
#include "CellData.h"
#include "massFlowRate.h"
#include "linearInterpolation.h"

Matrix::Matrix
(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const BoundaryConditions& boundaryConds
) : allFaces(faces), 
    allCells(cells),
    bcManager(boundaryConds)
{
    for (const auto& patch : bcManager.getPatches())
    {
        for (size_t i = patch.firstFaceIndex; i <= patch.lastFaceIndex; ++i)
        {
            faceToPatchMap[i] = &patch;
        }
    }
}

void Matrix::clear()
{
    tripletList.clear();

    size_t numCells = allCells.size();

    A_matrix.resize(numCells, numCells);
    A_matrix.setZero();

    b_vector.resize(numCells);
    b_vector.setZero();

    A_matrix.makeCompressed();
}

void Matrix::buildMatrix
(
    const ScalarField& phi,
    const ScalarField& phi_source,
    const VectorField& U_field,
    const ScalarField& Gamma,
    const ConvectionScheme& convScheme,
    const GradientScheme& gradScheme,
    const std::string& fieldName
)
{
    clear();

    size_t numCells = allCells.size();

    reserveTripletList();

    FaceFluxField mDotFace = 
        calculateMassFlowRate
        (
            allFaces,
            U_field,
            bcManager,
            std::map<size_t, const BoundaryPatch*>()
        );

    // Add source term
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        b_vector(cellIdx) += phi_source[cellIdx];
    }
    
    size_t numFaces = allFaces.size();
    for (size_t faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Face& face = allFaces[faceIdx];
        const size_t ownerIdx = face.ownerCell;
        
        Vector S_f = face.normal * face.area;
            
        if (face.isBoundary())
        {
            auto patchIt = faceToPatchMap.find(face.id);
            if (patchIt == faceToPatchMap.end()) 
            {
                std::cerr << "ERROR: Boundary face " << face.id 
                          << " not found in patch map! Owner cell " << ownerIdx
                          << " will have zero diagonal." << std::endl;
                continue;
            }
            
            const BoundaryData* bc = 
                bcManager.getFieldBC
                (
                    patchIt->second->patchName,
                    fieldName
                );
                  
            if (!bc) 
            {
                std::cerr << "ERROR: No boundary condition found for patch"
                          << " " << patchIt->second->patchName 
                          << " field " << fieldName << " face " << face.id
                          << " owner cell " << ownerIdx << std::endl;
                continue;
            }

            // Get geometric properties from owner perspective
            const Vector e_Pf = face.e_Pf;
            const Scalar d_Pf_mag = face.d_Pf_mag;

            Vector E_f = dot(S_f, e_Pf) * e_Pf;  //E_f = orthogonal component

            Scalar Gamma_f = Gamma[ownerIdx];

            Scalar a_diff = 
                Gamma_f * E_f.magnitude() / (d_Pf_mag + vSmallValue);

            if (bc->type == BCType::FIXED_VALUE)
            {
                Scalar phi_b = 0.0;
                
                // For momentum equations, extract component from vector BC
                if (fieldName == "U_x")
                {
                    phi_b = bc->vectorValue.x;
                }
                else if (fieldName == "U_y")
                {
                    phi_b = bc->vectorValue.y;
                }
                else if (fieldName == "U_z")
                {
                    phi_b = bc->vectorValue.z;
                }
                else  // For scalar fields (pressure, temperature, k, omega, etc.)
                {
                    phi_b = bc->scalarValue;
                }
                
                // Calculate mass flux based on boundary velocity
                Scalar a_conv = mDotFace[faceIdx];
                
                tripletList.emplace_back(ownerIdx, ownerIdx, a_diff);
                b_vector(ownerIdx) += a_diff * phi_b;
                b_vector(ownerIdx) -= a_conv * phi_b;

                // Non-orthogonal correction
                Vector T_f = S_f - E_f;

                Vector grad_phi_P = 
                    gradScheme.CellGradient(ownerIdx, phi, allCells);

                Scalar nonOrthogonalFlux = Gamma_f * dot(grad_phi_P, T_f);

                b_vector(ownerIdx) -= nonOrthogonalFlux;
            }
            else if (bc->type == BCType::ZERO_GRADIENT)
            {
                // Zero normal gradient -> no diffusive flux, only convection
                // Use cell center velocity for mass flux calculation
                Scalar a_conv = mDotFace[faceIdx];
                
                tripletList.emplace_back(ownerIdx, ownerIdx, a_conv);
            }
            else if (bc->type == BCType::NO_SLIP)
            {   // Simple Dirichlet BC: U = (0, 0, 0) at the wall
               
                tripletList.emplace_back(ownerIdx, ownerIdx, a_diff);

                // Non-orthogonal correction
                Vector T_f = S_f - E_f;

                Vector grad_phi_P = 
                    gradScheme.CellGradient(ownerIdx, phi, allCells);

                Scalar nonOrthogonalFlux = Gamma_f * dot(grad_phi_P, T_f);

                b_vector(ownerIdx) -= nonOrthogonalFlux;
            }
            else
            {
                // Undefined or unhandled BC type - default to zero gradient
                std::cerr << "Warning: Undefined boundary condition type for field " 
                          << fieldName << " on patch " 
                          << faceToPatchMap.at(face.id)->patchName 
                          << ". Applying zero gradient." << std::endl;
                
                // Calculate convective contribution for zero gradient
                Scalar mDot_face = mDotFace[faceIdx];
                Scalar a_conv = mDot_face;
                
                tripletList.emplace_back(ownerIdx, ownerIdx, a_conv);
            }
        }
        else // Internal face
        {
            const size_t neighborIdx = face.neighbourCell.value();

            Vector d_PN = allCells[neighborIdx].centroid - allCells[ownerIdx].centroid;
            Scalar d_PN_magnitude = d_PN.magnitude();
            Vector e_PN = d_PN / (d_PN_magnitude + vSmallValue);
            Vector E_f = dot(S_f, e_PN) * e_PN;

            // Harmonic interpolation for diffusion coefficient
            Scalar d_Pf = face.d_Pf_mag;
            Scalar d_Nf = face.d_Nf_mag.value();
            Scalar Gamma_f = d_PN_magnitude / ((d_Pf/(Gamma[ownerIdx] + vSmallValue)) + (d_Nf/(Gamma[neighborIdx] + vSmallValue)));
            Scalar a_diff = Gamma_f * E_f.magnitude() / (d_PN_magnitude + vSmallValue);

            // Non-orthogonal correction
            Vector T_f = S_f - E_f;
            Vector grad_phi_P = gradScheme.CellGradient(ownerIdx, phi, allCells);
            Vector grad_phi_N = gradScheme.CellGradient(neighborIdx, phi, allCells);

            Vector grad_phi_f = 
                gradScheme.FaceGradient
                (
                    face.id,
                    grad_phi_P,
                    grad_phi_N,
                    phi,
                    allCells,
                    allFaces,
                    bcManager,
                    fieldName
                );

            Scalar nonOrthogonalFlux = Gamma_f * dot(grad_phi_f, T_f);

            Scalar mDot_face = mDotFace[faceIdx];
            
            // Get convection coefficients
            Scalar a_P_conv, a_N_conv;
            convScheme.getFluxCoefficients(mDot_face, a_P_conv, a_N_conv);
            
            // Matrix coefficients for owner and neighbor cells
            // Owner cell contributions
            tripletList.emplace_back(ownerIdx, ownerIdx, a_diff + a_P_conv);
            tripletList.emplace_back(ownerIdx, neighborIdx, -a_diff + a_N_conv);
            
            // Neighbor cell contributions (flip convection coefficients)
            Scalar mass_flux_neighbor = -mDot_face;
            Scalar a_P_conv_N, a_N_conv_N;
            convScheme.getFluxCoefficients(mass_flux_neighbor, a_P_conv_N, a_N_conv_N);
            
            tripletList.emplace_back(neighborIdx, neighborIdx, a_diff + a_P_conv_N);
            tripletList.emplace_back(neighborIdx, ownerIdx, -a_diff + a_N_conv_N);
            
            // Source term contribution: FluxV_f = -μ_f (∇V)_f · T_f + ṁ_f (u_f_highRes - u_f_upwind)
            Scalar FluxV_f = -nonOrthogonalFlux;  // Non-orthogonal correction: -μ_f (∇V)_f · T_f
            
            // High-order correction: ṁ_f (u_f_highRes - u_f_upwind)
            if 
            (
                const CentralDifferenceScheme* cds = 
                    dynamic_cast<const CentralDifferenceScheme*>(&convScheme)
            )
            {
                FluxV_f += cds->calculateCentralDifferenceCorrection
                    (
                        face,
                        phi,
                        grad_phi_f,
                        mDot_face
                    );
            }
            else if 
            (
                const SecondOrderUpwindScheme* sous = 
                    dynamic_cast<const SecondOrderUpwindScheme*>(&convScheme)
            )
            {
                FluxV_f += sous->calculateSecondOrderCorrection
                    (
                        face,
                        phi,
                        grad_phi_P,
                        grad_phi_N,        
                        mDot_face
                    );
            }
            
            // Apply source term to owner and neighbor cells
            b_vector(ownerIdx) -= FluxV_f;
            b_vector(neighborIdx) += FluxV_f;
        }
    }

    A_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

void Matrix::buildPressureCorrectionMatrix(
    const FaceFluxField& RhieChowFlowRate,
    const FaceFluxField& DUf
)
{
    clear();

    reserveTripletList();

    // Build mass imbalance RHS by iterating over cells and their faces (respecting face signs)
    for (size_t cellIdx = 0; cellIdx < allCells.size(); ++cellIdx) 
    {
        const Cell& cell = allCells[cellIdx];
        
        // Mass imbalance contribution for this cell: -∑(sign_f * ṁ_f)
        Scalar massImbalance = 0.0;
        for (size_t i = 0; i < cell.faceIndices.size(); ++i)
        {
            const size_t faceIdx = cell.faceIndices[i];
            const int sign = cell.faceSigns[i];
            massImbalance += sign * RhieChowFlowRate[faceIdx];
        }
        
        // RHS = -mass_imbalance (negative because ∇·(D ∇p') = -∇·u*)
        b_vector[cellIdx] = -massImbalance;
    }

    // Build matrix coefficients by iterating over faces
    for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
    {
        const Face& face = allFaces[faceIdx];
        const size_t ownerIdx = face.ownerCell;

        if (face.isBoundary())
        {
            const std::string patchName = faceToPatchMap.at(face.id)->patchName;
            const BoundaryData* bc = bcManager.getFieldBC(patchName, "p");

            // Fixed-value pressure BC: pressure correction = 0
            if (bc && bc->type == BCType::FIXED_VALUE)
            {
                // Apply pCorr = 0 at fixed pressure boundaries
                tripletList.emplace_back(ownerIdx, ownerIdx, DUf[faceIdx]);
            }
        }
        else // Internal face
        {
            const size_t neighborIdx = face.neighbourCell.value();
            
            // Diagonal contributions
            tripletList.emplace_back(ownerIdx, ownerIdx, DUf[faceIdx]);
            tripletList.emplace_back(neighborIdx, neighborIdx, DUf[faceIdx]);
            
            // Off-diagonal contributions
            tripletList.emplace_back(ownerIdx, neighborIdx, -DUf[faceIdx]);
            tripletList.emplace_back(neighborIdx, ownerIdx, -DUf[faceIdx]);
        }
    }

    A_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

// Apply under-relaxation to the assembled linear system (Ax=b)
// Diagonal: a_c <- a_c / alpha
// RHS:      b   <- b + ((1 - alpha)/alpha) * a_c_original * phi_prev
void Matrix::relax(Scalar alpha, const ScalarField& phi_prev)
{
    if (alpha <= 0.0)
    {
        throw std::runtime_error("Matrix::relax: alpha must be positive");
    }

    const int n = static_cast<int>(A_matrix.rows());

    if (phi_prev.size() != static_cast<size_t>(n))
    {
        throw std::runtime_error
        (
            "Matrix::relax: phi_prev size mismatch with matrix size"
        );
    }

    // Cache original diagonal entries
    std::vector<Scalar> original_diag(static_cast<size_t>(n));

    for (int i = 0; i < n; ++i) 
    {
        original_diag[static_cast<size_t>(i)] = A_matrix.coeff(i, i);
    }

    // Scale diagonal
    for (int i = 0; i < n; ++i)
    {
        A_matrix.coeffRef(i, i) = original_diag[static_cast<size_t>(i)] / alpha;
    }

    // Update RHS
    const Scalar factor = (S(1.0) - alpha) / alpha;

    for (int i = 0; i < n; ++i)
    {
        b_vector(i) += factor * original_diag[static_cast<size_t>(i)] * phi_prev[static_cast<size_t>(i)];
    }
}

void Matrix::reserveTripletList()
{
    // Reserve for tripletList
    size_t internalFaces = 0, boundaryFaces = 0;

    for (const auto &f : allFaces) 
    {
        if (f.isBoundary()) ++boundaryFaces; else ++internalFaces;
    }
    
    size_t reserveSize = 4 * internalFaces + boundaryFaces;
    tripletList.reserve(reserveSize);
}