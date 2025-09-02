#include "massFlowRate.h"

FaceFluxField calculateMassFlowRate
(
    const std::vector<Face>& faces,
    const VectorField& U_field,
    const BoundaryConditions& bcManager,
    const std::map<size_t, const BoundaryPatch*>& /* faceToPatchMap */
) 
{
    FaceFluxField mdot("massFlowRate", faces.size(), 0.0);
    
    for (size_t faceIdx = 0; faceIdx < faces.size(); ++faceIdx)
    {
        const Face& face = faces[faceIdx];
        
        size_t ownerIdx = face.ownerCell;
        Vector S_f = face.normal * face.area;
        
        if (face.isBoundary())
        {
            // Use centralized boundary condition handling
            Vector U_face = 
                bcManager.calculateBoundaryFaceVectorValue(face, U_field, "U");
                
            mdot[faceIdx] = dot(U_face, S_f);
        }
        else
        {
            size_t neighborIdx = face.neighbourCell.value();
            
            // Distances from face centroid to cell centroids
            Scalar d_P = face.d_Pf_mag;
            Scalar d_N = face.d_Nf_mag.value();
            Scalar denom = d_P + d_N;

            // Weight for owner cell (note: inverted compared to distance)
            Scalar w_P = d_N / (denom + vSmallValue);
            
            // Interpolated velocity
            Vector U_f = w_P * U_field[ownerIdx] + (1.0 - w_P) * U_field[neighborIdx];
            
            mdot[faceIdx] = dot(U_f, S_f);
        }
    }
    
    return mdot;
}