#include "massFlowRate.h"

FaceFluxField calculateMassFlowRate(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const VectorField& U_field,
    Scalar rho
) {
    FaceFluxField mdot("massFlowRate", faces.size(), 0.0);
    
    for (size_t faceId = 0; faceId < faces.size(); ++faceId) {
        const Face& face = faces[faceId];
        
        // Skip faces without calculated geometric properties
        if (!face.geometricPropertiesCalculated) {
            mdot[faceId] = 0.0;
            continue;
        }
        
        size_t P = face.ownerCell;
        Vector S_f = face.normal * face.area;
        
        if (face.isBoundary()) {
            // Boundary face: use owner cell velocity only
            mdot[faceId] = rho * dot(U_field[P], S_f);
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