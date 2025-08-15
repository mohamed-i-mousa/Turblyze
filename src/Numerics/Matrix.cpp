#include "Matrix.h"
#include "KOmegaSST.h"
#include <iostream>
#include <algorithm>
#include "CellData.h"
#include "massFlowRate.h"

// ---- Private helpers ---- //
std::string Matrix::resolveBCFieldName(const std::string& fieldName) const {
    if (fieldName == "U_x" || fieldName == "U_y" || fieldName == "U_z") {
        return std::string("U");
    }
    if (fieldName == "p_prime") {
        return std::string("p");
    }
    return fieldName;
}

Scalar Matrix::computeDirichletValue(const BoundaryData* bc, const std::string& fieldName) const {
    if (!bc) return 0.0;
    const bool isUx = fieldName == "U_x";
    const bool isUy = fieldName == "U_y";
    const bool isUz = fieldName == "U_z";
    if (isUx || isUy || isUz) {
        if (bc->type == BCType::NO_SLIP) {
            return 0.0;
        }
        if (bc->type == BCType::FIXED_VALUE) {
            if (bc->valueType == BCValueType::VECTOR) {
                if (isUx) return bc->vectorValue.x;
                if (isUy) return bc->vectorValue.y;
                return bc->vectorValue.z;
            } else if (bc->valueType == BCValueType::SCALAR) {
                return bc->scalarValue;
            }
        }
        return 0.0;
    }
    // Scalar equations
    return bc->getFixedScalarValue();
}


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
    gradP_f = gradScheme.interpolateGradientsToFaces(gradP, p, allCells, allFaces, bcManager, "p");

    // Interpolate velocity gradients to faces
    gradUx_f = gradScheme.interpolateGradientsToFaces(gradUx, Ux, allCells, allFaces, bcManager, "U");
    gradUy_f = gradScheme.interpolateGradientsToFaces(gradUy, Uy, allCells, allFaces, bcManager, "U");
    gradUz_f = gradScheme.interpolateGradientsToFaces(gradUz, Uz, allCells, allFaces, bcManager, "U");

    // Compute turbulence gradients if turbulence model is available
    if (turbulenceModel) {
        const ScalarField& k_field = turbulenceModel->getK();
        const ScalarField& omega_field = turbulenceModel->getOmega();
        
        gradk = gradScheme.LeastSquares(k_field, allCells);
        gradOmega = gradScheme.LeastSquares(omega_field, allCells);
        
        // Interpolate turbulence gradients to faces
        gradk_f = gradScheme.interpolateGradientsToFaces(gradk, k_field, allCells, allFaces, bcManager, "k");
        gradOmega_f = gradScheme.interpolateGradientsToFaces(gradOmega, omega_field, allCells, allFaces, bcManager, "omega");
    }

    cachesValid = true;
}

void Matrix::buildMomentumMatrix(
    const std::string& fieldName,
    const ScalarField& phi,
    const ScalarField& phi_old,
    const ScalarField& phi_source,
    Scalar rho,
    const ScalarField& Gamma,
    TimeScheme timeScheme,
    Scalar dt,
    Scalar theta,
    const VectorField& grad_phi,
    const FaceVectorField& grad_phi_f,
    const ConvectionScheme& convScheme)
{
    clear();
    size_t numCells = allCells.size();

    // Use cached mass flow rate for momentum equations
    const FaceFluxField& mdot = mdotFaces;

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
                // const BoundaryData* bc = bcManager.getFieldBC(faceToPatchMap.at(face.id)->patchName, fieldName);
                const std::string bcField = resolveBCFieldName(fieldName);
                const BoundaryData* bc = bcManager.getFieldBC(faceToPatchMap.at(face.id)->patchName, bcField);
                if (!bc) continue;

                const Vector& e_Pf = face.e_Pf;
                Vector E_f = dot(S_f, e_Pf) * e_Pf;  // orthogonal component E_f (minimum approach)
                Scalar Gamma_f = Gamma[P];
                Scalar a_diff = Gamma_f * E_f.magnitude() / face.d_Pf_mag;
                Scalar a_P_conv, a_N_conv;
                convScheme.getFluxCoefficients(mdot[face.id], a_P_conv, a_N_conv);

                // if (bc->type == BCType::FIXED_VALUE) {
                if (bc->type == BCType::FIXED_VALUE) {
                    Scalar phi_b = computeDirichletValue(bc, fieldName);
                    // Include convective diagonal contribution for stability (TO BE CHECKED)
                    tripletList.emplace_back(P, P, a_diff + a_P_conv);
                    // RHS gets opposing side contribution (TO BE CHECKED)
                    b_vector(P) += (a_diff - a_N_conv) * phi_b;
                    // Non-orthogonal correction
                    Vector T_f = S_f - E_f;
                    Scalar flux_nonOrth = Gamma_f * dot(grad_phi[P], T_f);
                    b_vector(P) -= flux_nonOrth;
                } else if (bc->type == BCType::FIXED_GRADIENT) {
                    b_vector(P) -= Gamma_f * bc->getFixedScalarGradient() * E_f.magnitude();
                    tripletList.emplace_back(P, P, a_P_conv);
                    // Non-orthogonal correction
                    Vector T_f = S_f - E_f;
                    Scalar flux_nonOrth = Gamma_f * dot(grad_phi[P], T_f);
                    b_vector(P) -= flux_nonOrth;
                } else if (bc->type == BCType::ZERO_GRADIENT) {
                    // Zero normal gradient -> no diffusive flux; keep convective diag contribution only
                    tripletList.emplace_back(P, P, a_P_conv);
                    // No non-orthogonal correction for zero gradient
                }
            } 
            
            // ----- INTERNAL FACE LOGIC ----- //
            else { 
                size_t N = face.neighbourCell.value();
                Vector d_PN = allCells[N].centroid - allCells[P].centroid;
                Vector S_f = face.normal * face.area;
                Vector e_PN = d_PN / d_PN.magnitude();
                Vector E_f = dot(S_f, e_PN) * e_PN;  // orthogonal component
                // Harmonic interpolation of Gamma to the face
                Scalar d_P = (face.centroid - allCells[P].centroid).magnitude();
                Scalar d_N = (face.centroid - allCells[N].centroid).magnitude();
                Scalar d_PN_mag = d_PN.magnitude();
                Scalar denom = d_P / (Gamma[P] + 1e-20) + d_N / (Gamma[N] + 1e-20);
                Scalar Gamma_f = d_PN_mag / (denom + 1e-20);
                Scalar a_diff = Gamma_f * E_f.magnitude() / d_PN_mag;

                // Non-orthogonal correction
                Vector T_f = S_f - dot(S_f, e_PN) * e_PN; // perpendicular component of S_f
                Scalar flux_nonOrth = Gamma_f * dot(grad_phi_f[face.id], T_f);
                b_vector(P) -= flux_nonOrth;   // subtract from owner
                b_vector(N) += flux_nonOrth;   // add to neighbour (conservation)

                Scalar F = mdot[face.id];
                Scalar a_P_conv, a_N_conv;
                convScheme.getFluxCoefficients(F, a_P_conv, a_N_conv);

                // High-order correction
                Scalar highOrderCorrection = 0.0;
                if (const CentralDifferenceScheme* cds = dynamic_cast<const CentralDifferenceScheme*>(&convScheme)) {
                    highOrderCorrection = cds->calculateCentralDifferenceCorrection(face, phi, grad_phi_f, F);
                } else if (const SecondOrderUpwindScheme* sous = dynamic_cast<const SecondOrderUpwindScheme*>(&convScheme)) {
                    highOrderCorrection = sous->calculateSecondOrderCorrection(face, phi, grad_phi, F);
                } else if (dynamic_cast<const UpwindScheme*>(&convScheme)) {
                    highOrderCorrection = 0.0;
                }
                
                b_vector(P) -= highOrderCorrection;   // subtract from owner
                b_vector(N) += highOrderCorrection;   // add to neighbour (conservation)

                
                // Contribution to cell P
                tripletList.emplace_back(P, P, a_diff + a_P_conv);
                tripletList.emplace_back(P, N, -a_diff + a_N_conv);

                // Contribution to cell N is opposite
                tripletList.emplace_back(N, N, a_diff - a_N_conv);
                tripletList.emplace_back(N, P, -a_diff - a_P_conv);
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
                const Vector& e_Pf = face.e_Pf;
                Vector E_f = dot(S_f, e_Pf) * e_Pf;  // orthogonal component
                Scalar Gamma_f = Gamma[P];
                Scalar a_diff = Gamma_f * E_f.magnitude() / face.d_Pf_mag;

                // --- IMPLICIT PART (contributes to A and b) ---
                Scalar F_new = mdot[face.id];
                Scalar a_P_conv_new, a_N_conv_new; // a_N is not used for BCs but required by function
                convScheme.getFluxCoefficients(F_new, a_P_conv_new, a_N_conv_new);

                if (bc->type == BCType::FIXED_VALUE) {
                    // The implicit flux is: (D + a_P_conv_new)*phi_P - (D - min(F,0))*phi_b
                    // The phi_P term goes to the matrix diagonal
                    tripletList.emplace_back(P, P, theta * (a_diff + a_P_conv_new));
                    // The phi_b term goes to the source vector
                    b_vector(P) += theta * (a_diff - a_N_conv_new) * bc->getFixedScalarValue();
                } else if (bc->type == BCType::FIXED_GRADIENT) {
                    b_vector(P) -= theta * Gamma_f * bc->getFixedScalarGradient() * E_f.magnitude();
                    tripletList.emplace_back(P, P, theta * a_P_conv_new);
                }

                // --- EXPLICIT PART (contributes only to b) ---
                if (dt > 0 && theta < 1.0) {
                    Scalar F_old = mdot[face.id]; // Assuming U is constant over dt
                    Scalar a_P_conv_old, a_N_conv_old;
                    convScheme.getFluxCoefficients(F_old, a_P_conv_old, a_N_conv_old);

                    if (bc->type == BCType::FIXED_VALUE) {
                        // Total spatial flux from the previous time step
                        Scalar flux_old = (a_diff + a_P_conv_old) * phi_old[P] + (-a_diff + a_N_conv_old) * bc->getFixedScalarValue();
                        b_vector(P) -= (S(1.0) - theta) * flux_old;
                    } else if (bc->type == BCType::FIXED_GRADIENT) {
                    b_vector(P) -= (S(1.0) - theta) * Gamma_f * bc->getFixedScalarGradient() * E_f.magnitude();
                        b_vector(P) += (S(1.0) - theta) * (a_P_conv_old * phi_old[P]);
                    }
                }

                // ----- Non-orthogonal correction (over-relaxed) ----- //
                Vector t_f_b = S_f - dot(S_f, e_Pf) * e_Pf;
                Scalar flux_nonOrth_bnd = Gamma_f * dot(grad_phi[P], t_f_b);
                b_vector(P) -= flux_nonOrth_bnd;
            } else {
                // --- INTERNAL FACE LOGIC ---
                size_t N = face.neighbourCell.value();
                Vector d_PN = allCells[N].centroid - allCells[P].centroid;
                Vector S_f = face.normal * face.area;
                Vector e_PN = d_PN / d_PN.magnitude();
                Vector E_f = dot(S_f, e_PN) * e_PN;  // orthogonal component
                Scalar d_P = (face.centroid - allCells[P].centroid).magnitude();
                Scalar d_N = (face.centroid - allCells[N].centroid).magnitude();
                Scalar d_PN_mag = d_PN.magnitude();
                Scalar denom = d_P / (Gamma[P] + 1e-20) + d_N / (Gamma[N] + 1e-20);
                Scalar Gamma_f = d_PN_mag / (denom + 1e-20);
                Scalar a_diff = Gamma_f * E_f.magnitude() / d_PN_mag;

                // Implicit part using pre-computed mass flow rate
                Scalar F_new = mdot[face.id];
                Scalar a_P_conv_new, a_N_conv_new;
                convScheme.getFluxCoefficients(F_new, a_P_conv_new, a_N_conv_new);
                
                // Contribution to P
                tripletList.emplace_back(P, P, theta * (a_diff + a_P_conv_new));
                tripletList.emplace_back(P, N, theta * (-a_diff + a_N_conv_new));
                
                // Contribution to N is opposite
                tripletList.emplace_back(N, N, theta * (a_diff - a_N_conv_new));
                tripletList.emplace_back(N, P, theta * (-a_diff - a_P_conv_new));

                // Explicit part
                if (dt > 0 && theta < 1.0) {
                    Scalar F_old = mdot[face.id]; // Assuming U is constant over dt
                    Scalar a_P_conv_old, a_N_conv_old;
                    convScheme.getFluxCoefficients(F_old, a_P_conv_old, a_N_conv_old);
                    
                    Scalar flux_old = (a_diff + a_P_conv_old) * phi_old[P] + (-a_diff + a_N_conv_old) * phi_old[N];

                    b_vector(P) -= (S(1.0) - theta) * flux_old;
                    b_vector(N) += (S(1.0) - theta) * flux_old;
                }
            }
        }
    }

    A_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

void Matrix::buildScalarTransportMatrix(
    const std::string& fieldName,
    const ScalarField& phi,
    const ScalarField& phi_old,
    const VectorField& U_field,
    const ScalarField& phi_source,
    Scalar rho,
    const ScalarField& Gamma,
    TimeScheme timeScheme,
    Scalar dt,
    Scalar theta,
    const ConvectionScheme& convScheme)
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
    else {
        // For other fields, compute gradient locally
        grad_phi = gradScheme.LeastSquares(phi, allCells);
        grad_phi_ptr = &grad_phi;
    }
    
    const VectorField& grad_phi_ref = *grad_phi_ptr;

    // Use cached face gradients for known fields, otherwise compute locally
    FaceVectorField grad_phi_f("grad_phi_f", allFaces.size(), Vector(0.0,0.0,0.0));
    if      (isUx)          grad_phi_f = gradUx_f;
    else if (isUy)          grad_phi_f = gradUy_f;
    else if (isUz)          grad_phi_f = gradUz_f;
    else if (isPressure)    grad_phi_f = gradP_f;
    else if (isk)           grad_phi_f = gradk_f;
    else if (isOmega)       grad_phi_f = gradOmega_f;
    else {
        // For other fields, interpolate gradients to faces
        grad_phi_f = gradScheme.interpolateGradientsToFaces(grad_phi, phi, allCells, allFaces, bcManager, fieldName);
    }

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
                const std::string bcField = resolveBCFieldName(fieldName);
                const BoundaryData* bc = bcManager.getFieldBC(faceToPatchMap.at(face.id)->patchName, bcField);
                if (!bc) continue;

                const Vector& e_Pf = face.e_Pf;
                Vector E_f = dot(S_f, e_Pf) * e_Pf;  // orthogonal component E_f (minimum approach)
                Scalar Gamma_f = Gamma[P];
                Scalar a_diff = Gamma_f * E_f.magnitude() / face.d_Pf_mag;
                Scalar a_P_conv, a_N_conv;
                convScheme.getFluxCoefficients(mdot[face.id], a_P_conv, a_N_conv);

                if (bc->type == BCType::FIXED_VALUE) {
                    Scalar phi_b = computeDirichletValue(bc, fieldName);
                    tripletList.emplace_back(P, P, a_diff);
                    b_vector(P) += (a_diff - mdot[face.id]) * phi_b;
                    // Non-orthogonal correction
                    Vector T_f = S_f - E_f;
                    Scalar flux_nonOrth = Gamma_f * dot(grad_phi_ref[P], T_f);
                    b_vector(P) -= flux_nonOrth;
                } else if (bc->type == BCType::FIXED_GRADIENT) {
                    b_vector(P) -= Gamma_f * bc->getFixedScalarGradient() * E_f.magnitude();
                    tripletList.emplace_back(P, P, a_P_conv);
                    // Non-orthogonal correction
                    Vector T_f = S_f - E_f;
                    Scalar flux_nonOrth = Gamma_f * dot(grad_phi_ref[P], T_f);
                    b_vector(P) -= flux_nonOrth;
                } else if (bc->type == BCType::ZERO_GRADIENT) {
                    // Zero normal gradient -> no diffusive flux; keep convective diag contribution only
                    tripletList.emplace_back(P, P, a_P_conv);
                    // No non-orthogonal correction for zero gradient
                }
            } 
            
            // ----- INTERNAL FACE LOGIC ----- //
            else { 
                size_t N = face.neighbourCell.value();
                Vector d_PN = allCells[N].centroid - allCells[P].centroid;
                Vector S_f = face.normal * face.area;
                Vector e_PN = d_PN / d_PN.magnitude();
                Vector E_f = dot(S_f, e_PN) * e_PN;  // orthogonal component
                Scalar d_P = (face.centroid - allCells[P].centroid).magnitude();
                Scalar d_N = (face.centroid - allCells[N].centroid).magnitude();
                // Harmonic interpolation of Gamma at face
                Scalar denom_gamma = d_P / (Gamma[P] + 1e-20) + d_N / (Gamma[N] + 1e-20);
                Scalar Gamma_f = (d_P + d_N) / (denom_gamma + 1e-20);
                Scalar a_diff = Gamma_f * E_f.magnitude() / d_PN.magnitude();

                // Non-orthogonal correction
                Vector T_f = S_f - dot(S_f, e_PN) * e_PN; // perpendicular component of S_f
                Scalar flux_nonOrth = Gamma_f * dot(grad_phi_f[face.id], T_f);
                b_vector(P) -= flux_nonOrth;   // subtract from owner
                b_vector(N) += flux_nonOrth;   // add to neighbour (conservation)

                Scalar F = mdot[face.id];
                Scalar a_P_conv, a_N_conv;
                convScheme.getFluxCoefficients(F, a_P_conv, a_N_conv);

                // High-order correction
                Scalar highOrderCorrection = 0.0;
                if (const CentralDifferenceScheme* cds = dynamic_cast<const CentralDifferenceScheme*>(&convScheme)) {
                    highOrderCorrection = cds->calculateCentralDifferenceCorrection(face, phi, grad_phi_f, F);
                } else if (const SecondOrderUpwindScheme* sous = dynamic_cast<const SecondOrderUpwindScheme*>(&convScheme)) {
                    highOrderCorrection = sous->calculateSecondOrderCorrection(face, phi, grad_phi_ref, F);
                } else if (dynamic_cast<const UpwindScheme*>(&convScheme)) {
                    highOrderCorrection = 0.0;
                }
                
                b_vector(P) -= highOrderCorrection;   // subtract from owner
                b_vector(N) += highOrderCorrection;   // add to neighbour (conservation)

                
                // Contribution to cell P
                tripletList.emplace_back(P, P, a_diff + a_P_conv);
                tripletList.emplace_back(P, N, -a_diff + a_N_conv);

                // Contribution to cell N is opposite
                tripletList.emplace_back(N, N, a_diff - a_N_conv);
                tripletList.emplace_back(N, P, -a_diff - a_P_conv);
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
                const Vector& e_Pf = face.e_Pf;
                Vector E_f = dot(S_f, e_Pf) * e_Pf;  // orthogonal component
                Scalar Gamma_f = Gamma[P];
                Scalar a_diff = Gamma_f * E_f.magnitude() / face.d_Pf_mag;

                // --- IMPLICIT PART (contributes to A and b) ---
                Scalar F_new = mdot[face.id];
                Scalar a_P_conv_new, a_N_conv_new; // a_N is not used for BCs but required by function
                convScheme.getFluxCoefficients(F_new, a_P_conv_new, a_N_conv_new);

                if (bc->type == BCType::FIXED_VALUE) {
                    // The implicit flux is: (D + a_P_conv_new)*phi_P - (D - min(F,0))*phi_b
                    // The phi_P term goes to the matrix diagonal
                    tripletList.emplace_back(P, P, theta * (a_diff + a_P_conv_new));
                    // The phi_b term goes to the source vector
                    b_vector(P) += theta * (a_diff - a_N_conv_new) * bc->getFixedScalarValue();
                } else if (bc->type == BCType::FIXED_GRADIENT) {
                    b_vector(P) -= theta * Gamma_f * bc->getFixedScalarGradient() * E_f.magnitude();
                    tripletList.emplace_back(P, P, theta * a_P_conv_new);
                }

                // --- EXPLICIT PART (contributes only to b) ---
                if (dt > 0 && theta < 1.0) {
                    Scalar F_old = mdot[face.id]; // Assuming U is constant over dt
                    Scalar a_P_conv_old, a_N_conv_old;
                    convScheme.getFluxCoefficients(F_old, a_P_conv_old, a_N_conv_old);

                    if (bc->type == BCType::FIXED_VALUE) {
                        // Total spatial flux from the previous time step
                        Scalar flux_old = (a_diff + a_P_conv_old) * phi_old[P] + (-a_diff + a_N_conv_old) * bc->getFixedScalarValue();
                        b_vector(P) -= (S(1.0) - theta) * flux_old;
                    } else if (bc->type == BCType::FIXED_GRADIENT) {
                    b_vector(P) -= (S(1.0) - theta) * Gamma_f * bc->getFixedScalarGradient() * E_f.magnitude();
                        b_vector(P) += (S(1.0) - theta) * (a_P_conv_old * phi_old[P]);
                    }
                }

                // ----- Non-orthogonal correction (over-relaxed) ----- //
                Vector t_f_b = S_f - dot(S_f, e_Pf) * e_Pf;
                Scalar flux_nonOrth_bnd = Gamma_f * dot(grad_phi_ref[P], t_f_b);
                b_vector(P) -= flux_nonOrth_bnd;
            } else {
                // --- INTERNAL FACE LOGIC ---
                size_t N = face.neighbourCell.value();
                Vector d_PN = allCells[N].centroid - allCells[P].centroid;
                Vector S_f = face.normal * face.area;
                Vector e_PN = d_PN / d_PN.magnitude();
                Vector E_f = dot(S_f, e_PN) * e_PN;  // orthogonal component
                Scalar d_P = (face.centroid - allCells[P].centroid).magnitude();
                Scalar d_N = (face.centroid - allCells[N].centroid).magnitude();
                // Harmonic interpolation of Gamma at face
                Scalar denom_gamma = d_P / (Gamma[P] + 1e-20) + d_N / (Gamma[N] + 1e-20);
                Scalar Gamma_f = (d_P + d_N) / (denom_gamma + 1e-20);
                Scalar a_diff = Gamma_f * E_f.magnitude() / d_PN.magnitude();

                // Implicit part using pre-computed mass flow rate
                Scalar a_P_conv_new, a_N_conv_new;
                convScheme.getFluxCoefficients(mdot[face.id], a_P_conv_new, a_N_conv_new);
                
                // Contribution to P
                tripletList.emplace_back(P, P, theta * (a_diff + a_P_conv_new));
                tripletList.emplace_back(P, N, theta * (-a_diff + a_N_conv_new));
                
                // Contribution to N is opposite
                tripletList.emplace_back(N, N, theta * (a_diff - a_N_conv_new));
                tripletList.emplace_back(N, P, theta * (-a_diff - a_P_conv_new));

                // Explicit part
                if (dt > 0 && theta < 1.0) {
                    Scalar F_old = mdot[face.id]; // Assuming U is constant over dt
                    Scalar a_P_conv_old, a_N_conv_old;
                    convScheme.getFluxCoefficients(F_old, a_P_conv_old, a_N_conv_old);
                    
                    Scalar flux_old = (a_diff + a_P_conv_old) * phi_old[P] + (-a_diff + a_N_conv_old) * phi_old[N];

                    b_vector(P) -= (S(1.0) - theta) * flux_old;
                    b_vector(N) += (S(1.0) - theta) * flux_old;
                }
            }
        }
    }

    A_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

void Matrix::buildPressureMatrix(
    const FaceFluxField& massFlux,
    const ScalarField& a_Ux,
    const ScalarField& a_Uy,
    const ScalarField& a_Uz,
    Scalar rho)
{
    clear();
    const size_t numCells = allCells.size();

    // Build RHS from mass imbalance: b = -sum(massFlux)
    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx) {
        const Cell& cell = allCells[cellIdx];
        Scalar mass_imbalance = 0.0;
        for (size_t j = 0; j < cell.faceIndices.size(); ++j) {
            const size_t faceIdx = cell.faceIndices[j];
            const int sign = cell.faceSigns[j];
            mass_imbalance += sign * massFlux[faceIdx];
        }
        b_vector(cellIdx) = -mass_imbalance;
    }

    // Detect whether there is a fixed pressure boundary
    bool hasFixedPressureBC = false;
    for (const auto& patch : bcManager.patches) {
        const BoundaryData* bc = bcManager.getFieldBC(patch.patchName, "p");
        if (bc && bc->type == BCType::FIXED_VALUE) { hasFixedPressureBC = true; break; }
    }

    // Reserve triplets (rough estimate)
    size_t internalFaces = 0, boundaryFaces = 0;
    for (const auto &f : allFaces) {
        if (f.isBoundary()) ++boundaryFaces; else ++internalFaces;
    }
    tripletList.clear();
    tripletList.reserve(4 * internalFaces + boundaryFaces);

    // Assemble diffusion-like matrix using Rhie-Chow-consistent coefficients
    for (const auto& face : allFaces) {
        if (!face.geometricPropertiesCalculated) continue;

        const size_t P = face.ownerCell;
        if (face.isBoundary()) {
            // Fixed pressure boundary: enforce p' = 0
            const std::string patchName = faceToPatchMap.at(face.id)->patchName;
            const BoundaryData* bc = bcManager.getFieldBC(patchName, "p");
            if (bc && bc->type == BCType::FIXED_VALUE) {
                tripletList.emplace_back(P, P, 1e10);
                b_vector(P) = 0.0;
            }
            continue;
        }

        const size_t N = face.neighbourCell.value();
        const Vector d_PN = allCells[N].centroid - allCells[P].centroid;
        const Scalar d_P = face.d_Pf_mag;
        const Scalar d_N = face.d_Nf_mag.value();
        const Scalar total_dist = d_P + d_N;
        const Scalar w_P = d_N / total_dist;
        const Scalar w_N = d_P / total_dist;

        // Per-component diagonal at face
        const Scalar a_face_x = w_P * a_Ux[P] + w_N * a_Ux[N];
        const Scalar a_face_y = w_P * a_Uy[P] + w_N * a_Uy[N];
        const Scalar a_face_z = w_P * a_Uz[P] + w_N * a_Uz[N];

        // Minimum-correction coefficient with non-orthogonality (pressure equation needs rho)
        const Vector S_f = face.normal * face.area;
        const Scalar denom = d_PN.magnitudeSquared();
        const Scalar numer = d_PN.x * S_f.x / (a_face_x + 1e-20)
                           + d_PN.y * S_f.y / (a_face_y + 1e-20)
                           + d_PN.z * S_f.z / (a_face_z + 1e-20);
        const Scalar D_f = rho * (numer / denom);

        tripletList.emplace_back(P, P, D_f);
        tripletList.emplace_back(P, N, -D_f);
        tripletList.emplace_back(N, N, D_f);
        tripletList.emplace_back(N, P, -D_f);
    }

    A_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

    // Anchor if no fixed-pressure boundary is present
    if (!hasFixedPressureBC && A_matrix.rows() > 0) {
        A_matrix.coeffRef(0, 0) += 1e12;
        b_vector(0) = 0.0;
    }
}

// Apply under-relaxation to the assembled linear system (Ax=b)
// Diagonal: a_c <- a_c / alpha
// RHS:      b   <- b + ((1 - alpha)/alpha) * a_c_original * phi_prev
void Matrix::relax(Scalar alpha, const ScalarField& phi_prev)
{
    if (alpha <= 0.0) {
        throw std::runtime_error("Matrix::relax: alpha must be positive");
    }

    const int n = static_cast<int>(A_matrix.rows());
    if (phi_prev.size() != static_cast<size_t>(n)) {
        throw std::runtime_error("Matrix::relax: phi_prev size mismatch with matrix size");
    }

    // Cache original diagonal entries
    std::vector<Scalar> original_diag(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        original_diag[static_cast<size_t>(i)] = A_matrix.coeff(i, i);
    }

    // Scale diagonal
    for (int i = 0; i < n; ++i) {
        A_matrix.coeffRef(i, i) = original_diag[static_cast<size_t>(i)] / alpha;
    }

    // Update RHS
    const Scalar factor = (S(1.0) - alpha) / alpha;
    for (int i = 0; i < n; ++i) {
        b_vector(i) += factor * original_diag[static_cast<size_t>(i)] * phi_prev[static_cast<size_t>(i)];
    }
}