#include "massFlowRate.h"

FaceFluxField calculateMassFlowRate(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const VectorField& U_field,
    Scalar rho,
    const BoundaryConditions& bcManager,
    const std::map<size_t, const BoundaryPatch*>& faceToPatchMap
) {
    FaceFluxField mdot("massFlowRate", faces.size(), 0.0);
    
    for (size_t faceId = 0; faceId < faces.size(); ++faceId) {
        const Face& face = faces[faceId];
        
        size_t P = face.ownerCell;
        Vector S_f = face.normal * face.area;
        
        if (face.isBoundary()) {
            Vector U_face;
            
            // Get boundary condition information
            const BoundaryPatch* patch = faceToPatchMap.at(faceId);
            const BoundaryData* bc = bcManager.getFieldBC(patch->patchName, "U");
            
            // Check if boundary condition is properly configured
            if (bc == nullptr) {
                // No boundary condition configured, use zero gradient as fallback
                std::cerr << "Warning: No boundary condition configured for patch '" 
                          << patch->patchName << "' field 'U'. Using zero gradient fallback." << std::endl;
                U_face = U_field[P];
            } else {
                switch (bc->type) {
                case BCType::FIXED_VALUE:
                    U_face = bc->vectorValue;
                    break;
                    
                case BCType::FIXED_GRADIENT: {
                    Vector d_Pf = face.centroid - cells[P].centroid;
                    U_face = U_field[P] + bc->vectorGradient * d_Pf.magnitude();
                    break;
                }
                    
                case BCType::ZERO_GRADIENT:
                    U_face = U_field[P];
                    break;
                    
                case BCType::NO_SLIP:
                    U_face = Vector(0.0, 0.0, 0.0);
                    break;
                    
                default:
                    // Fallback to owner cell velocity
                    U_face = U_field[P];
                    break;
                }
            }
            
            mdot[faceId] = rho * dot(U_face, S_f);
        } else {
            // Internal face: distance-weighted interpolation
            size_t N = face.neighbourCell.value();
            
            // Calculate distances from face centroid to cell centroids
            Scalar d_P = (face.centroid - cells[P].centroid).magnitude();
            Scalar d_N = (face.centroid - cells[N].centroid).magnitude();
            Scalar denom = d_P + d_N;
            Scalar w = d_P / denom;
            
            // Interpolated velocity: U_field[P] * (1 - w) + U_field[N] * w
            Vector U_f = U_field[P] * (1.0 - w) + U_field[N] * w;
            
            mdot[faceId] = rho * dot(U_f, S_f);
        }
    }
    
    return mdot;
}