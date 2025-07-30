#include "GradientScheme.h"

#include <stdexcept>
#include <eigen3/Eigen/Dense>


    
/* Cell-centered Least-Squares Gradient Reconstruction
* 
* For each cell P, we solve the overdetermined system:
*   phi_N - phi_P = ∇φ_P · (r_N - r_P) + O(|r_N - r_P|²)
* 
* This leads to the normal equations:  A^T A ∇φ_P = A^T b
* where:
*   A_ij = w_i * (r_i - r_P)_j
*   b_i = w_i * (phi_i - phi_P)
*   w_i = 1/|r_i - r_P|² (inverse distance squared weighting)
* 
* The solution ∇φ_P is obtained by solving the 3x3 system using LLT decomposition.
*/
VectorField GradientScheme::LeastSquares(
    const ScalarField& phi,
    const std::vector<Cell>& allCells) const
{
    if (allCells.empty()) {
        throw std::invalid_argument("LeastSquares: Empty cell list provided");
    }
    
    size_t numCells = allCells.size();
    if (phi.size() != numCells) {
        throw std::invalid_argument("LeastSquares: Field size (" + std::to_string(phi.size()) + 
                                   ") does not match cell count (" + std::to_string(numCells) + ")");
    }
    
    VectorField grad_phi("grad(" + phi.name + ")", numCells, Vector(0,0,0));

    // Pre-allocate Eigen matrices to avoid repeated allocations
    Eigen::Matrix<Scalar,3,3> ATA;
    Eigen::Matrix<Scalar,3,1> ATb;
    Eigen::Matrix<Scalar,3,1> r_vector;
    
    // Performance monitoring
    size_t cellsProcessed = 0;
    size_t cellsSkipped = 0;

    // Loop over all cells and solve small 3x3 least-squares system
    for (size_t i = 0; i < numCells; ++i) {
        const Cell& cell = allCells[i];

        if (cell.neighbourCellIndices.empty()) {
            cellsSkipped++;
            continue; 
        }

        if (cell.neighbourCellIndices.size() < 3) {
            cellsSkipped++;
            continue;
        }

        ATA.setZero();
        ATb.setZero();

        Scalar totalWeight = S(0.0);
        size_t validNeighbors = 0;

        // Process internal cell neighbors only
        for (size_t neighbourId : cell.neighbourCellIndices) {
            if (neighbourId >= numCells) continue; 
            
            const Cell& neighbour = allCells[neighbourId];
            Vector r = neighbour.centroid - cell.centroid;
            
            // Check for degenerate cases (cells too close together)
            Scalar r_mag_sq = r.magnitudeSquared();
            if (r_mag_sq < GRADIENT_TOLERANCE) continue;
            
            Scalar w = S(1.0) / r_mag_sq; // weight ∝ 1/|r|²
            totalWeight += w;
            validNeighbors++;

            r_vector << r.x, r.y, r.z;
            ATA.noalias() += w * (r_vector * r_vector.transpose());
            
            Scalar delta_phi = phi[neighbourId] - phi[i];
            ATb.noalias() += w * delta_phi * r_vector;
        }

        if (validNeighbors < 3) {
            cellsSkipped++;
            continue;
        }
        
        cellsProcessed++;

        // Add small regularization to improve conditioning
        Scalar regularization = totalWeight * 1e-12;
        ATA(0,0) += regularization;
        ATA(1,1) += regularization;
        ATA(2,2) += regularization;

        // Efficient solver for 3x3 systems
        Eigen::LLT<Eigen::Matrix<Scalar,3,3>> llt(ATA);
        if (llt.info() == Eigen::Success) {
            Eigen::Matrix<Scalar,3,1> g = llt.solve(ATb);
            grad_phi[i] = Vector(g(0), g(1), g(2));
        } else {
            // Fallback to more robust solver if LLT fails
            Eigen::FullPivLU<Eigen::Matrix<Scalar,3,3>> lu(ATA);
            if (lu.isInvertible()) {
                Eigen::Matrix<Scalar,3,1> g = lu.solve(ATb);
                grad_phi[i] = Vector(g(0), g(1), g(2));
            }
        }
    }
    
    // Print statistics
    std::cout << "LeastSquares gradient: " << cellsProcessed << " cells processed, " 
              << cellsSkipped << " cells skipped" << std::endl;

    return grad_phi;
}

/* Interpolate cell-centered gradients to face values using the scheme:
* ∇φ_f = average(∇φ_f) + [(φ_N - φ_P)/2 - dot(average(∇φ_f), e_PN)] * e_PN
* 
* where:
*   average(∇φ_f) = (g_P * ∇φ_P + g_N * ∇φ_N)
*   e_PN = d_PN / |d_PN|
*   g_P, g_N are interpolation weights (distance-based)
* 
* This scheme provides second-order accuracy and ensures consistency
* between the gradient interpolation and the underlying field values.
*/

FaceVectorField GradientScheme::interpolateGradientsToFaces(
    const VectorField& grad_phi,
    const ScalarField& phi,
    const std::vector<Cell>& allCells,
    const std::vector<Face>& allFaces) const
{
    // Input validation
    if (allFaces.empty()) {
        throw std::invalid_argument("interpolateGradientsToFaces: Empty face list provided");
    }
    
    if (allCells.empty()) {
        throw std::invalid_argument("interpolateGradientsToFaces: Empty cell list provided");
    }
    
    size_t numFaces = allFaces.size();
    size_t numCells = allCells.size();
    
    if (grad_phi.size() != numCells) {
        throw std::invalid_argument("interpolateGradientsToFaces: Gradient field size (" + 
                                   std::to_string(grad_phi.size()) + 
                                   ") does not match cell count (" + std::to_string(numCells) + ")");
    }
    
    if (phi.size() != numCells) {
        throw std::invalid_argument("interpolateGradientsToFaces: Scalar field size (" + 
                                   std::to_string(phi.size()) + 
                                   ") does not match cell count (" + std::to_string(numCells) + ")");
    }
    
    // Initialize face gradient field
    FaceVectorField grad_phi_faces("grad_phi_faces", numFaces, Vector(0,0,0));
    
    // Performance monitoring
    size_t internalFacesProcessed = 0;
    size_t boundaryFacesProcessed = 0;
    size_t facesSkipped = 0;
    
    // Loop over all faces
    for (size_t faceId = 0; faceId < numFaces; ++faceId) {
        const Face& face = allFaces[faceId];
        
        if (!face.geometricPropertiesCalculated) {
            facesSkipped++;
            continue;
        }
        
        size_t P = face.ownerCell;
        
        if (face.isBoundary()) {
            // For boundary faces, use the owner cell gradient
            // This is a simple approach; more sophisticated boundary treatments
            // could be implemented based on boundary conditions
            grad_phi_faces[faceId] = grad_phi[P];
            boundaryFacesProcessed++;
        } else {
            // For internal faces, use the specified interpolation scheme
            size_t N = face.neighbourCell.value();
            
            // Calculate the vector from P to N
            Vector d_PN = allCells[N].centroid - allCells[P].centroid;
            Scalar d_PN_mag = d_PN.magnitude();
            
            if (d_PN_mag < GRADIENT_TOLERANCE) {
                // Cells are too close, use simple average
                grad_phi_faces[faceId] = S(0.5) * (grad_phi[P] + grad_phi[N]);
                continue;
            }
            
            Vector e_PN = d_PN / d_PN_mag;
            
            // Calculate interpolation weights (distance-based)
            Scalar d_Pf = (face.centroid - allCells[P].centroid).magnitude();
            Scalar d_Nf = (face.centroid - allCells[N].centroid).magnitude();
            Scalar total_dist = d_Pf + d_Nf;
            
            Scalar g_P, g_N;
            if (total_dist < GRADIENT_TOLERANCE) {
                // Face is equidistant, use simple average
                g_P = S(0.5);
                g_N = S(0.5);
            } else {
                // Distance-weighted interpolation
                g_P = d_Nf / total_dist; // Weight for cell P
                g_N = d_Pf / total_dist; // Weight for cell N
            }
            
            // Calculate average gradient at face
            Vector grad_avg = g_P * grad_phi[P] + g_N * grad_phi[N];
            
            // Calculate the correction term
            Scalar phi_diff = phi[N] - phi[P];
            Scalar correction = (phi_diff / d_PN_mag) - dot(grad_avg, e_PN);
            
            // Apply the interpolation scheme
            grad_phi_faces[faceId] = grad_avg + correction * e_PN;
            internalFacesProcessed++;
        }
    }
    
    // Print performance statistics
    std::cout << "Gradient interpolation: " << internalFacesProcessed << " internal faces, " 
              << boundaryFacesProcessed << " boundary faces, " << facesSkipped << " faces skipped" << std::endl;
    
    return grad_phi_faces;
}