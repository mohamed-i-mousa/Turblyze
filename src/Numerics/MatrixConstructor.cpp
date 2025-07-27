#include "MatrixConstructor.h"
#include <iostream>
#include <algorithm>
#include "CellData.h"

MatrixConstructor::MatrixConstructor(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const BoundaryConditions& boundaryConds,
    const GradientScheme& gradientScheme
) : allFaces(faces), allCells(cells), bcManager(boundaryConds), gradScheme(gradientScheme)
{
    // Build the face-to-patch map once for efficient boundary lookup
    for (const auto& patch : bcManager.patches) {
        if(patch.lastFaceIndex >= patch.firstFaceIndex) {
            for (size_t i = patch.firstFaceIndex; i <= patch.lastFaceIndex; ++i) {
                faceToPatchMap[i] = &patch;
            }
        }
    }
}

void MatrixConstructor::clear() {
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

void MatrixConstructor::constructScalarTransportMatrix(
    const std::string& fieldName,
    const ScalarField& phi,
    const ScalarField& phi_old,
    const VectorField& U_field,
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

    // Pre-compute gradient field
    VectorField grad_phi = gradScheme.LeastSquares(phi, allCells);

    // Rough estimate to minimise reallocations of tripletList
    size_t internalFaces = 0, boundaryFaces = 0;
    for (const auto &f : allFaces) {
        if (f.isBoundary()) ++boundaryFaces; else ++internalFaces;
    }
    size_t reserveSize = 4 * internalFaces + 2 * boundaryFaces + numCells; // diag/time terms
    tripletList.reserve(reserveSize);

    // ----- STEADY-STATE a*phi = b ----- //
    if (timeScheme == TimeScheme::Steady) {
        for (const auto& face : allFaces) {
            if (!face.geometricPropertiesCalculated) continue;

            size_t P = face.ownerCell;
            Vector S_f = face.normal * face.area;

            if (face.isBoundary()) {
                const BoundaryData* bc = bcManager.getFieldBC(faceToPatchMap.at(face.id)->patchName, fieldName);
                if (!bc) continue;

                Scalar D = Gamma * face.area / (face.centroid - allCells[P].centroid).magnitude();
                Scalar F = rho * dot(U_field[P], S_f);
                Scalar a_P_conv, a_N_conv;
                convScheme.getFluxCoefficients(F, a_P_conv, a_N_conv, face, allCells);

                if (bc->type == BCType::FIXED_VALUE) {
                    tripletList.emplace_back(P, P, D);
                    b_vector(P) += (D - F) * bc->getFixedScalarValue();
                } else if (bc->type == BCType::FIXED_GRADIENT) {                        // TODO: check if this is correct
                    b_vector(P) -= Gamma * bc->getFixedScalarGradient() * face.area;
                    tripletList.emplace_back(P, P, a_P_conv);
                }

                // ----- Non-orthogonal correction (over-relaxed) ----- //
                Vector d_Pf = face.centroid - allCells[P].centroid;
                Vector e_Pf = d_Pf.normalized();
                Vector T_f = S_f - dot(S_f, e_Pf) * e_Pf;
                Scalar flux_nonOrth = Gamma * dot(grad_phi[P], T_f);
                b_vector(P) -= flux_nonOrth;
            } else { // Internal Face
                size_t N = face.neighbourCell.value();
                Vector d_PN = allCells[N].centroid - allCells[P].centroid;
                Vector S_f = face.normal * face.area;
                Scalar D = Gamma * face.area / d_PN.magnitude();

                // Distance-weighted interpolation of the face velocity
                Scalar d_P = (face.centroid - allCells[P].centroid).magnitude();
                Scalar d_N = (face.centroid - allCells[N].centroid).magnitude();
                Scalar denom = d_P + d_N;
                Scalar w_P = d_N / denom;
                Vector U_f = w_P * U_field[P] + (1.0 - w_P) * U_field[N];

                // Non-orthogonal correction
                Vector e_PN = d_PN / d_PN.magnitude();
                Vector T_f = S_f - dot(S_f, e_PN) * e_PN; // perpendicular component of S_f
                Vector grad_f = w_P * grad_phi[P] + (1.0 - w_P) * grad_phi[N];
                Scalar flux_nonOrth = Gamma * dot(grad_f, T_f);
                b_vector(P) -= flux_nonOrth;   // subtract from owner
                b_vector(N) += flux_nonOrth;   // add to neighbour (conservation)

                Scalar F = rho * dot(U_f, S_f);
                Scalar a_P_conv, a_N_conv;
                convScheme.getFluxCoefficients(F, a_P_conv, a_N_conv, face, allCells);

                // Contribution to cell P
                tripletList.emplace_back(P, P, D + a_P_conv);
                tripletList.emplace_back(P, N, -D + a_N_conv);

                // Contribution to cell N is opposite
                tripletList.emplace_back(N, N, D - a_N_conv);
                tripletList.emplace_back(N, P, -D - a_P_conv);
            }
        }
    } else { // TRANSIENT
            // Add the base time derivative term
    if (dt > 0) {
        for (size_t i = 0; i < numCells; ++i) {
            Scalar term = rho * allCells[i].volume / dt;
            tripletList.emplace_back(i, i, term);
            b_vector(i) += term * phi_old[i];
        }
    }

    // Loop over all faces to add spatial terms
    for (const auto& face : allFaces) {
        if (!face.geometricPropertiesCalculated) continue;

        size_t P = face.ownerCell;

        if (face.isBoundary()) {
            // --- BOUNDARY FACE LOGIC ---
            const BoundaryData* bc = bcManager.getFieldBC(faceToPatchMap.at(face.id)->patchName, fieldName);
            if (!bc) continue; // Skip if no BC is set for this field

            Vector S_f = face.normal * face.area;
            Scalar D = Gamma * face.area / (face.centroid - allCells[P].centroid).magnitude();

            // --- IMPLICIT PART (contributes to A and b) ---
            Scalar F_new = rho * dot(U_field[P], S_f);
            Scalar a_P_conv_new, a_N_conv_new; // a_N is not used for BCs but required by function
            convScheme.getFluxCoefficients(F_new, a_P_conv_new, a_N_conv_new, face, allCells);

            if (bc->type == BCType::FIXED_VALUE) {
                // The implicit flux is: (D + a_P_conv_new)*phi_P - (D - min(F,0))*phi_b
                // The phi_P term goes to the matrix diagonal
                tripletList.emplace_back(P, P, theta * (D + a_P_conv_new));
                // The phi_b term goes to the source vector
                b_vector(P) += theta * (D - a_N_conv_new) * bc->getFixedScalarValue();
            } else if (bc->type == BCType::FIXED_GRADIENT) {
                b_vector(P) -= theta * Gamma * bc->getFixedScalarGradient() * face.area;
                tripletList.emplace_back(P, P, theta * a_P_conv_new);
            }

            // --- EXPLICIT PART (contributes only to b) ---
            if (dt > 0 && theta < 1.0) {
                Scalar F_old = rho * dot(U_field[P], S_f); // Assuming U is constant over dt
                Scalar a_P_conv_old, a_N_conv_old;
                convScheme.getFluxCoefficients(F_old, a_P_conv_old, a_N_conv_old, face, allCells);

                if (bc->type == BCType::FIXED_VALUE) {
                    // Total spatial flux from the previous time step
                    Scalar flux_old = (D + a_P_conv_old) * phi_old[P] + (-D + a_N_conv_old) * bc->getFixedScalarValue();
                    b_vector(P) -= (1.0 - theta) * flux_old;
                } else if (bc->type == BCType::FIXED_GRADIENT) {
                     b_vector(P) -= (1.0 - theta) * Gamma * bc->getFixedScalarGradient() * face.area;
                     b_vector(P) += (1.0 - theta) * (a_P_conv_old * phi_old[P]);
                }
            }

            // ----- Non-orthogonal correction (over-relaxed) ----- //
            Vector d_Pf = face.centroid - allCells[P].centroid;
            Vector e_Pf = d_Pf.normalized();
            Vector t_f_b = S_f - dot(S_f, e_Pf) * e_Pf;
            Scalar flux_nonOrth_bnd = Gamma * dot(grad_phi[P], t_f_b);
            b_vector(P) -= flux_nonOrth_bnd;
        } else {
            // --- INTERNAL FACE LOGIC ---
            size_t N = face.neighbourCell.value();
            Vector d_PN = allCells[N].centroid - allCells[P].centroid;
            Vector S_f = face.normal * face.area;
            Scalar D = Gamma * face.area / d_PN.magnitude();

            // Implicit part – distance-weighted face velocity
            Scalar d_P = (face.centroid - allCells[P].centroid).magnitude();
            Scalar d_N = (face.centroid - allCells[N].centroid).magnitude();
            Scalar denom = d_P + d_N;
            Scalar w_P = d_N / denom;
            Vector U_f_new = w_P * U_field[P] + (1.0 - w_P) * U_field[N];
            Scalar F_new = rho * dot(U_f_new, S_f);
            Scalar a_P_conv_new, a_N_conv_new;
            convScheme.getFluxCoefficients(F_new, a_P_conv_new, a_N_conv_new, face, allCells);
            
            // Contribution to P
            tripletList.emplace_back(P, P, theta * (D + a_P_conv_new));
            tripletList.emplace_back(P, N, theta * (-D + a_N_conv_new));
            
            // Contribution to N is opposite
            tripletList.emplace_back(N, N, theta * (D - a_N_conv_new));
            tripletList.emplace_back(N, P, theta * (-D - a_P_conv_new));

            // Explicit part
            if (dt > 0 && theta < 1.0) {
                Vector U_f_old = U_f_new; // Assuming U is constant over dt
                Scalar F_old = rho * dot(U_f_old, S_f);
                Scalar a_P_conv_old, a_N_conv_old;
                convScheme.getFluxCoefficients(F_old, a_P_conv_old, a_N_conv_old, face, allCells);
                
                Scalar flux_old = (D + a_P_conv_old) * phi_old[P] + (-D + a_N_conv_old) * phi_old[N];

                b_vector(P) -= (1.0 - theta) * flux_old;
                b_vector(N) += (1.0 - theta) * flux_old;
            }
        }
    }
    }

    A_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}