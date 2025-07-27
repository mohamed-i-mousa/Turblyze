#include "GradientScheme.h"

#include <stdexcept>
#include <eigen3/Eigen/Dense>

    
Scalar GradientScheme::getBoundaryFaceValue(
    const Face& face,
    size_t ownerCellId,
    const ScalarField& phi,
    const std::vector<Cell>& allCells,
    const BoundaryData& bc) const
{
    switch (bc.type) {
        case BCType::FIXED_VALUE:
            // Directly return the fixed value for the scalar field
            if (bc.valueType == BCValueType::SCALAR) {
                return bc.scalarValue;
            }
            // Add logic for vector fields if needed, e.g., for U, V, W components
            break;
            
        case BCType::FIXED_GRADIENT: {
            // phi_f = phi_P + (grad_phi . n) * |d_Pf|
            // This is an approximation. We assume face normal is parallel to d_Pf.
            const Vector d_Pf = face.centroid - allCells[ownerCellId].centroid;
            return phi[ownerCellId] + bc.scalarGradient * d_Pf.magnitude();
        }

        case BCType::ZERO_GRADIENT:
        
        case BCType::NO_SLIP: // For scalar fields like pressure, no-slip on wall is often zero-gradient
        
        default:
            return phi[ownerCellId];
    }
    // Fallback for unhandled cases
    return phi[ownerCellId];
}

VectorField GradientScheme::GreenGauss(
    const ScalarField& phi,
    const std::vector<Face>& allFaces,
    const std::vector<Cell>& allCells,
    const BoundaryConditions& bcManager) const
{
    size_t numCells = allCells.size();
    VectorField grad_phi("grad(" + phi.name + ")", numCells);

    // Create a temporary map from face ID to its boundary patch for efficient lookup
    std::map<size_t, const BoundaryPatch*> faceToPatchMap;
    for (const auto& patch : bcManager.patches) {
        for (size_t i = patch.firstFaceIndex; i <= patch.lastFaceIndex; ++i) {
            faceToPatchMap[i] = &patch;
        }
    }

    // First, compute phi values on all faces
    FaceData<Scalar> phi_f("phi_faces", allFaces.size());
    for(size_t i = 0; i < allFaces.size(); ++i) {
        const auto& face = allFaces[i];
        if (face.isBoundary()) {
            auto it = faceToPatchMap.find(face.id);
            if (it != faceToPatchMap.end()) {
                const BoundaryPatch* patch = it->second;
                const BoundaryData* bc = bcManager.getFieldBC(patch->patchName, phi.name);
                if (bc) {
                    phi_f[i] = getBoundaryFaceValue(face, face.ownerCell, phi, allCells, *bc);
                } else { // Fallback: zero-gradient
                    phi_f[i] = phi[face.ownerCell];
                }
            }
        } else {
            // For internal faces, use linear interpolation (central differencing)
            const Vector& P = allCells[face.ownerCell].centroid;
            const Vector& N = allCells[face.neighbourCell.value()].centroid;
            Scalar w = (face.centroid - P).magnitude() / (N - P).magnitude();
            phi_f[i] = phi[face.ownerCell] * (S(1.0) - w) + phi[face.neighbourCell.value()] * w;
        }
    }
    
    // Now, calculate gradients for each cell
    for (size_t i = 0; i < numCells; ++i) {
        const Cell& cell = allCells[i];
        Vector sum_phi_f_Sf(S(0.0), S(0.0), S(0.0));

        for (size_t j = 0; j < cell.faceIndices.size(); ++j) {
            size_t face_id = cell.faceIndices[j];
            const Face& face = allFaces[face_id];
            
            // Use face area vector with correct orientation for the cell
            Vector S_f = face.normal * face.area * S(cell.faceSigns[j]);
            
            sum_phi_f_Sf += phi_f[face_id] * S_f;
        }
        
        if (cell.volume > 1e-12) {
            grad_phi[i] = sum_phi_f_Sf / cell.volume;
        }
    }
    return grad_phi;
}

// -----------------------------------------------------------------------------
// Least Squares gradient reconstruction
// -----------------------------------------------------------------------------
VectorField GradientScheme::LeastSquares(
    const ScalarField& phi,
    const std::vector<Cell>& allCells) const
{
    size_t numCells = allCells.size();
    VectorField grad_phi("grad(" + phi.name + ")", numCells, Vector(0,0,0));

    if (numCells == 0) return grad_phi;

    // Loop over all cells and solve small 3x3 least-squares system
    for (size_t i = 0; i < numCells; ++i) {
        const Cell& cell = allCells[i];

        if (cell.neighbourCellIndices.empty()) {
            continue; // leave zero gradient if no neighbours available
        }

        Eigen::Matrix<Scalar,3,3> ATA = Eigen::Matrix<Scalar,3,3>::Zero();
        Eigen::Matrix<Scalar,3,1> ATb = Eigen::Matrix<Scalar,3,1>::Zero();

        for (size_t neighbourId : cell.neighbourCellIndices) {
            if (neighbourId >= numCells) continue; // safety
            const Cell& neighbour = allCells[neighbourId];
            Vector r = neighbour.centroid - cell.centroid;
            Scalar w = S(1.0) / (r.magnitudeSquared()); // weight ∝ 1/|r|²

            Eigen::Matrix<Scalar,3,1> r_vector;
            r_vector << r.x, r.y, r.z;

            ATA += w * (r_vector * r_vector.transpose());            // ATA is A^T * A
            Scalar delta_phi = phi[neighbourId] - phi[i];
            ATb += w * delta_phi * r_vector;                         // ATb is A^T * b
        }

        // Solve ATA * g = ATb; LDLT for robustness
        Eigen::FullPivLU<Eigen::Matrix<Scalar,3,3>> lu(ATA);
        if (lu.isInvertible()) {
            Eigen::Matrix<Scalar,3,1> g = lu.solve(ATb);
            grad_phi[i] = Vector(g(0), g(1), g(2));
        }
    }

    return grad_phi;
}