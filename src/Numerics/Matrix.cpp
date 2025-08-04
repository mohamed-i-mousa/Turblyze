#include "Matrix.h"
#include "KOmegaSST.h"
#include <iostream>
#include <algorithm>
#include "CellData.h"
#include "massFlowRate.h"

Matrix::Matrix(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const BoundaryConditions& boundaryConds,
    const GradientScheme& gradientScheme
) : allFaces(faces), allCells(cells), bcManager(boundaryConds), gradScheme(gradientScheme),
      gradP("gradP", cells.size(), Vector(0.0,0.0,0.0)),
      gradUx("gradUx", cells.size(), Vector(0.0,0.0,0.0)),
      gradUy("gradUy", cells.size(), Vector(0.0,0.0,0.0)),
      gradUz("gradUz", cells.size(), Vector(0.0,0.0,0.0)),
      gradk("gradk", cells.size(), Vector(0.0,0.0,0.0)),
      gradOmega("gradOmega", cells.size(), Vector(0.0,0.0,0.0)),
      mdotFaces("mdot", faces.size(), 0.0),
      gradP_f("gradP_f", faces.size(), Vector(0.0,0.0,0.0)),
      gradUx_f("gradUx_f", faces.size(), Vector(0.0,0.0,0.0)),
      gradUy_f("gradUy_f", faces.size(), Vector(0.0,0.0,0.0)),
      gradUz_f("gradUz_f", faces.size(), Vector(0.0,0.0,0.0)),
      gradk_f("gradk_f", faces.size(), Vector(0.0,0.0,0.0)),
      gradOmega_f("gradOmega_f", faces.size(), Vector(0.0,0.0,0.0))
{
    // Build the face-to-patch map once for efficient boundary lookup
    for (const auto& patch : bcManager.patches) {
        for (size_t i = patch.firstFaceIndex; i <= patch.lastFaceIndex; ++i) {
            faceToPatchMap[i] = &patch;
        }
    }
}

void Matrix::clear() {
    tripletList.clear();
    size_t numCells = allCells.size();
    if (numCells > 0) {
        A_matrix.resize(numCells, numCells);
        A_matrix.setZero();
        b_vector.resize(numCells);
        b_vector.setZero();
    } else {
        A_matrix.resize(0, 0);
        b_vector.resize(0);
    }
    A_matrix.makeCompressed();
}

void Matrix::refreshIterationCaches(const PressureField& p, const VelocityField& U, Scalar rho, const KOmegaSST* turbulenceModel)
{
    // Compute cell-centre pressure gradient once
    gradP = gradScheme.LeastSquares(p, allCells);

    // Compute velocity component gradients once
    // Extract scalar fields for LeastSquares helper
    ScalarField Ux("Ux", U.size());
    ScalarField Uy("Uy", U.size());
    ScalarField Uz("Uz", U.size());
    for (size_t c = 0; c < U.size(); ++c) {
        Ux[c] = U[c].x;
        Uy[c] = U[c].y;
        Uz[c] = U[c].z;
    }
    gradUx = gradScheme.LeastSquares(Ux, allCells);
    gradUy = gradScheme.LeastSquares(Uy, allCells);
    gradUz = gradScheme.LeastSquares(Uz, allCells);

    // Face-level mdot
    mdotFaces = calculateMassFlowRate(allFaces, allCells, U, rho, bcManager, faceToPatchMap);

    // Interpolate pressure gradient to faces
    gradP_f = gradScheme.interpolateGradientsToFaces(gradP, p, allCells, allFaces);

    // Interpolate velocity gradients to faces
    gradUx_f = gradScheme.interpolateGradientsToFaces(gradUx, Ux, allCells, allFaces);
    gradUy_f = gradScheme.interpolateGradientsToFaces(gradUy, Uy, allCells, allFaces);
    gradUz_f = gradScheme.interpolateGradientsToFaces(gradUz, Uz, allCells, allFaces);

    // Compute turbulence gradients if turbulence model is available
    if (turbulenceModel) {
        const ScalarField& k_field = turbulenceModel->getK();
        const ScalarField& omega_field = turbulenceModel->getOmega();
        
        gradk = gradScheme.LeastSquares(k_field, allCells);
        gradOmega = gradScheme.LeastSquares(omega_field, allCells);
        
        // Interpolate turbulence gradients to faces
        gradk_f = gradScheme.interpolateGradientsToFaces(gradk, k_field, allCells, allFaces);
        gradOmega_f = gradScheme.interpolateGradientsToFaces(gradOmega, omega_field, allCells, allFaces);
    }

    cachesValid = true;
}

void Matrix::constructScalarTransportMatrix(
    const std::string& fieldName,
    const ScalarField& phi,
    const ScalarField& phi_old,
    const VectorField& U_field,
    const ScalarField& phi_source,
    Scalar rho,
    Scalar Gamma,
    TimeScheme timeScheme,
    Scalar dt,
    Scalar theta,
    const ConvectionDiscretization& convScheme)
{
    clear();
    size_t numCells = allCells.size();
    if (numCells == 0) return;

    VectorField grad_phi("grad_tmp", allCells.size(), Vector(0.0,0.0,0.0));
    const bool isUx = fieldName == "U_x";
    const bool isUy = fieldName == "U_y";
    const bool isUz = fieldName == "U_z";
    const bool isPressure = fieldName == "p" || fieldName == "p_prime";
    const bool isk = fieldName == "k";
    const bool isOmega = fieldName == "omega";

    // Use a pointer to avoid unnecessary copies of gradient fields
    const VectorField* grad_phi_ptr = nullptr;
    if      (isUx)          grad_phi_ptr = &gradUx;
    else if (isUy)          grad_phi_ptr = &gradUy;
    else if (isUz)          grad_phi_ptr = &gradUz;
    else if (isPressure)    grad_phi_ptr = &gradP;
    else if (isk)           grad_phi_ptr = &gradk;
    else if (isOmega)       grad_phi_ptr = &gradOmega;

    
    const VectorField& grad_phi_ref = *grad_phi_ptr;

    // Use cached face gradients for known fields, otherwise compute locally
    FaceVectorField grad_phi_f("grad_phi_f", allFaces.size(), Vector(0.0,0.0,0.0));
    if      (isUx)          grad_phi_f = gradUx_f;
    else if (isUy)          grad_phi_f = gradUy_f;
    else if (isUz)          grad_phi_f = gradUz_f;
    else if (isPressure)    grad_phi_f = gradP_f;
    else if (isk)           grad_phi_f = gradk_f;
    else if (isOmega)       grad_phi_f = gradOmega_f;

    // Re-use cached mdot for velocity equations; otherwise compute locally
    FaceFluxField mdot = mdotFaces; // copy cached by default
    if (!(isUx || isUy || isUz)) {
        mdot = calculateMassFlowRate(allFaces, allCells, U_field, rho, bcManager, faceToPatchMap);
    }

    // Rough estimate to minimise reallocations of tripletList
    size_t internalFaces = 0, boundaryFaces = 0;
    for (const auto &f : allFaces) {
        if (f.isBoundary()) ++boundaryFaces; else ++internalFaces;
    }
    size_t reserveSize = 4 * internalFaces + 2 * boundaryFaces + numCells; // diag/time terms
    tripletList.reserve(reserveSize);

    // ----- STEADY-STATE ----- //
    if (timeScheme == TimeScheme::Steady) {
        // Add source term
        for (size_t i = 0; i < numCells; ++i) {
            b_vector(i) += phi_source[i];
        }
        
        for (const auto& face : allFaces) {
            if (!face.geometricPropertiesCalculated) continue;

            size_t P = face.ownerCell;
            Vector S_f = face.normal * face.area;

            // ----- BOUNDARY FACE LOGIC ----- //
            if (face.isBoundary()) {
                const BoundaryData* bc = bcManager.getFieldBC(faceToPatchMap.at(face.id)->patchName, fieldName);
                if (!bc) continue;

                Vector d_Pf = face.centroid - allCells[P].centroid;
                Vector e_Pf = d_Pf / d_Pf.magnitude();
                Vector E_f = dot(S_f, e_Pf) * e_Pf;  // orthogonal component E_f (minimum approach)
                Scalar D = Gamma * E_f.magnitude() / d_Pf.magnitude();
                Scalar F = mdot[face.id];
                Scalar a_P_conv, a_N_conv;
                convScheme.getFluxCoefficients(F, a_P_conv, a_N_conv);

                if (bc->type == BCType::FIXED_VALUE) {
                    tripletList.emplace_back(P, P, D);
                    b_vector(P) += (D - F) * bc->getFixedScalarValue();
                } else if (bc->type == BCType::FIXED_GRADIENT) {
                    b_vector(P) -= Gamma * bc->getFixedScalarGradient() * E_f.magnitude();
                    tripletList.emplace_back(P, P, a_P_conv);
                }

                // Non-orthogonal correction
                Vector T_f = S_f - E_f;
                Scalar flux_nonOrth = Gamma * dot(grad_phi_ref[P], T_f);
                b_vector(P) -= flux_nonOrth;
            } 
            
            // ----- INTERNAL FACE LOGIC ----- //
            else { 
                size_t N = face.neighbourCell.value();
                Vector d_PN = allCells[N].centroid - allCells[P].centroid;
                Vector S_f = face.normal * face.area;
                Vector e_PN = d_PN / d_PN.magnitude();
                Vector E_f = dot(S_f, e_PN) * e_PN;  // orthogonal component
                Scalar D = Gamma * E_f.magnitude() / d_PN.magnitude();

                // Non-orthogonal correction
                Vector T_f = S_f - dot(S_f, e_PN) * e_PN; // perpendicular component of S_f
                Scalar flux_nonOrth = Gamma * dot(grad_phi_f[face.id], T_f);
                b_vector(P) -= flux_nonOrth;   // subtract from owner
                b_vector(N) += flux_nonOrth;   // add to neighbour (conservation)

                Scalar F = mdot[face.id];
                Scalar a_P_conv, a_N_conv;
                convScheme.getFluxCoefficients(F, a_P_conv, a_N_conv);

                // High-order correction
                Scalar correction = 0.0;
                if (const CentralDifferenceScheme* cds = dynamic_cast<const CentralDifferenceScheme*>(&convScheme)) {
                    correction = cds->calculateCentralDifferenceCorrection(face, allCells, phi, grad_phi_f, F);
                } else if (const SecondOrderUpwindScheme* sous = dynamic_cast<const SecondOrderUpwindScheme*>(&convScheme)) {
                    correction = sous->calculateSecondOrderCorrection(face, allCells, phi, grad_phi_ref, F);
                } else if (dynamic_cast<const UpwindScheme*>(&convScheme)) {
                    correction = 0.0;
                }
                
                b_vector(P) -= correction;   // subtract from owner
                b_vector(N) += correction;   // add to neighbour (conservation)

                
                // Contribution to cell P
                tripletList.emplace_back(P, P, D + a_P_conv);
                tripletList.emplace_back(P, N, -D + a_N_conv);

                // Contribution to cell N is opposite
                tripletList.emplace_back(N, N, D - a_N_conv);
                tripletList.emplace_back(N, P, -D - a_P_conv);
            }
        }
    } else { // TRANSIENT

        // Add the time derivative term and pressure gradient term
        if (dt > 0) {
            for (size_t i = 0; i < numCells; ++i) {
                Scalar term = rho * allCells[i].volume / dt;    // implicit transient term
                tripletList.emplace_back(i, i, term);   
                b_vector(i) += term * phi_old[i];       // explicit transient term
                b_vector(i) += phi_source[i];           // source term 
            }
        }

        // Add the pressure gradient term
        // Loop over all faces to add spatial terms
        for (const auto& face : allFaces) {
            if (!face.geometricPropertiesCalculated) continue;

            size_t P = face.ownerCell;

            if (face.isBoundary()) {
                // --- BOUNDARY FACE LOGIC ---
                const BoundaryData* bc = bcManager.getFieldBC(faceToPatchMap.at(face.id)->patchName, fieldName);
                if (!bc) continue; // Skip if no BC is set for this field

                Vector S_f = face.normal * face.area;
                Vector d_Pf = face.centroid - allCells[P].centroid;
                Vector e_Pf = d_Pf / d_Pf.magnitude();
                Vector E_f = dot(S_f, e_Pf) * e_Pf;  // orthogonal component
                Scalar D = Gamma * E_f.magnitude() / d_Pf.magnitude();

                // --- IMPLICIT PART (contributes to A and b) ---
                Scalar F_new = mdot[face.id];
                Scalar a_P_conv_new, a_N_conv_new; // a_N is not used for BCs but required by function
                convScheme.getFluxCoefficients(F_new, a_P_conv_new, a_N_conv_new);

                if (bc->type == BCType::FIXED_VALUE) {
                    // The implicit flux is: (D + a_P_conv_new)*phi_P - (D - min(F,0))*phi_b
                    // The phi_P term goes to the matrix diagonal
                    tripletList.emplace_back(P, P, theta * (D + a_P_conv_new));
                    // The phi_b term goes to the source vector
                    b_vector(P) += theta * (D - a_N_conv_new) * bc->getFixedScalarValue();
                } else if (bc->type == BCType::FIXED_GRADIENT) {
                    b_vector(P) -= theta * Gamma * bc->getFixedScalarGradient() * E_f.magnitude();
                    tripletList.emplace_back(P, P, theta * a_P_conv_new);
                }

                // --- EXPLICIT PART (contributes only to b) ---
                if (dt > 0 && theta < 1.0) {
                    Scalar F_old = mdot[face.id]; // Assuming U is constant over dt
                    Scalar a_P_conv_old, a_N_conv_old;
                    convScheme.getFluxCoefficients(F_old, a_P_conv_old, a_N_conv_old);

                    if (bc->type == BCType::FIXED_VALUE) {
                        // Total spatial flux from the previous time step
                        Scalar flux_old = (D + a_P_conv_old) * phi_old[P] + (-D + a_N_conv_old) * bc->getFixedScalarValue();
                        b_vector(P) -= (S(1.0) - theta) * flux_old;
                    } else if (bc->type == BCType::FIXED_GRADIENT) {
                        b_vector(P) -= (S(1.0) - theta) * Gamma * bc->getFixedScalarGradient() * E_f.magnitude();
                        b_vector(P) += (S(1.0) - theta) * (a_P_conv_old * phi_old[P]);
                    }
                }

                // ----- Non-orthogonal correction (over-relaxed) ----- //
                Vector t_f_b = S_f - dot(S_f, e_Pf) * e_Pf;
                Scalar flux_nonOrth_bnd = Gamma * dot(grad_phi_ref[P], t_f_b);
                b_vector(P) -= flux_nonOrth_bnd;
            } else {
                // --- INTERNAL FACE LOGIC ---
                size_t N = face.neighbourCell.value();
                Vector d_PN = allCells[N].centroid - allCells[P].centroid;
                Vector S_f = face.normal * face.area;
                Vector e_PN = d_PN / d_PN.magnitude();
                Vector E_f = dot(S_f, e_PN) * e_PN;  // orthogonal component
                Scalar D = Gamma * E_f.magnitude() / d_PN.magnitude();

                // Implicit part using pre-computed mass flow rate
                Scalar F_new = mdot[face.id];
                Scalar a_P_conv_new, a_N_conv_new;
                convScheme.getFluxCoefficients(F_new, a_P_conv_new, a_N_conv_new);
                
                // Contribution to P
                tripletList.emplace_back(P, P, theta * (D + a_P_conv_new));
                tripletList.emplace_back(P, N, theta * (-D + a_N_conv_new));
                
                // Contribution to N is opposite
                tripletList.emplace_back(N, N, theta * (D - a_N_conv_new));
                tripletList.emplace_back(N, P, theta * (-D - a_P_conv_new));

                // Explicit part
                if (dt > 0 && theta < 1.0) {
                    Scalar F_old = mdot[face.id]; // Assuming U is constant over dt
                    Scalar a_P_conv_old, a_N_conv_old;
                    convScheme.getFluxCoefficients(F_old, a_P_conv_old, a_N_conv_old);
                    
                    Scalar flux_old = (D + a_P_conv_old) * phi_old[P] + (-D + a_N_conv_old) * phi_old[N];

                    b_vector(P) -= (S(1.0) - theta) * flux_old;
                    b_vector(N) += (S(1.0) - theta) * flux_old;
                }
            }
        }
    }

    A_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}