#include "GradientScheme.h"

#include <stdexcept>
#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>


VectorField GradientScheme::LeastSquares
(
    const ScalarField& phi,
    const std::vector<Cell>& allCells
) const
{
    size_t numCells = allCells.size();

    if (phi.size() != numCells) 
    {
        throw std::invalid_argument
        (
            "LeastSquares: Field size (" + std::to_string(phi.size()) 
          + ") does not match cell count (" + std::to_string(numCells) + ")"
        );
    }
    
    VectorField grad_phi("grad(" + phi.name + ")", numCells, Vector(0,0,0));

    Eigen::Matrix<Scalar,3,3> ATA;
    Eigen::Matrix<Scalar,3,1> ATb;
    Eigen::Matrix<Scalar,3,1> r_vector;
    
    size_t cellsProcessed = 0;
    size_t cellsSkipped = 0;

    for (size_t i = 0; i < numCells; ++i) 
    {
        const Cell& cell = allCells[i];

        if (cell.neighbourCellIndices.size() < 3)
        {
            cellsSkipped++;
            continue;
        }

        ATA.setZero();
        ATb.setZero();

        Scalar totalWeight = S(0.0);
        size_t validNeighbors = 0;

        // Process internal cell neighbors only
        for (size_t neighbourId : cell.neighbourCellIndices) 
        {
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

        if (validNeighbors < 3) 
        {
            cellsSkipped++;
            continue;
        }
        
        cellsProcessed++;

        // Add small regularization to improve conditioning
        Scalar regularization = totalWeight * 1e-12;
        ATA(0,0) += regularization;
        ATA(1,1) += regularization;
        ATA(2,2) += regularization;

        Eigen::LLT<Eigen::Matrix<Scalar,3,3>> llt(ATA);
        if (llt.info() == Eigen::Success) 
        {
            Eigen::Matrix<Scalar,3,1> g = llt.solve(ATb);
            grad_phi[i] = Vector(g(0), g(1), g(2));
        }
        else
        {
            // Fallback to more robust solver if LLT fails
            Eigen::FullPivLU<Eigen::Matrix<Scalar,3,3>> lu(ATA);

            if (lu.isInvertible()) 
            {
                Eigen::Matrix<Scalar,3,1> g = lu.solve(ATb);
                grad_phi[i] = Vector(g(0), g(1), g(2));
            }
        }
    }

    return grad_phi;
}

FaceVectorField GradientScheme::interpolateGradientsToFaces
(
    const VectorField& grad_phi,
    const ScalarField& phi,
    const std::vector<Cell>& allCells,
    const std::vector<Face>& allFaces,
    const BoundaryConditions& boundaryConditions,
    const std::string& fieldName
) const
{
    
    size_t numFaces = allFaces.size();
    
    // Initialize face gradient field
    FaceVectorField grad_phi_faces("grad_phi_faces", numFaces, Vector(0,0,0));
    
    // Performance monitoring
    size_t internalFacesProcessed = 0;
    size_t boundaryFacesProcessed = 0;
    
    // Loop over all faces
    for (size_t faceId = 0; faceId < numFaces; ++faceId) 
    {
        const Face& face = allFaces[faceId];
        
        size_t P = face.ownerCell;
        
        if (face.isBoundary())
        {
            grad_phi_faces[faceId] = 
                calculateBoundaryFaceGradient
                (
                    face, grad_phi[P], phi, boundaryConditions, fieldName
                );

            boundaryFacesProcessed++;
        }
        else
        {
            // For internal faces, use the specified interpolation scheme
            size_t N = face.neighbourCell.value();
            
            // Calculate the vector from P to N
            Vector d_PN = allCells[N].centroid - allCells[P].centroid;
            Scalar d_PN_mag = d_PN.magnitude();
            
            if (d_PN_mag < GRADIENT_TOLERANCE)
            {
                // Cells are too close, use simple average
                grad_phi_faces[faceId] = S(0.5) * (grad_phi[P] + grad_phi[N]);
                continue;
            }
            
            Vector e_PN = d_PN / d_PN_mag;
            
            // Calculate interpolation weights (distance-based)
            Scalar d_Pf = face.d_Pf_mag;
            Scalar d_Nf = face.d_Nf_mag.value();
            Scalar total_dist = d_Pf + d_Nf;
            
            Scalar g_P, g_N;

            if (total_dist < GRADIENT_TOLERANCE)
            {
                throw std::runtime_error
                (
                    "GradientScheme::interpolateGradientsToFaces: "
                    "Face is equidistant, use simple average"
                );
            } 
            else 
            {
                // Distance-weighted interpolation
                g_P = d_Nf / total_dist;
                g_N = d_Pf / total_dist;
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

    return grad_phi_faces;
}

Vector GradientScheme::calculateBoundaryFaceGradient
(
    const Face& face,
    const Vector& cellGradient,
    const ScalarField& phi,
    const BoundaryConditions& boundaryConditions,
    const std::string& fieldName
) const
{
    // Find the boundary patch for this face
    const BoundaryPatch* patch = nullptr;
    for (const auto& p : boundaryConditions.patches)
    {
        if (face.id >= p.firstFaceIndex && face.id <= p.lastFaceIndex)
        {
            patch = &p;
            break;
        }
    }

    // Get the boundary condition for this field 
    const BoundaryData* bc = 
        boundaryConditions.getFieldBC(patch->patchName, fieldName);
    
    // If no BC specified, use cell gradient
    if (!bc) return cellGradient;
    
    switch (bc->type) {
        case BCType::FIXED_VALUE: {
            // For fixed value BC, calculate normal gradient from boundary-cell
            // difference and keep tangential components from cell gradient
            
            Scalar boundaryValue = S(0.0);
            if (bc->valueType == BCValueType::SCALAR) {
                boundaryValue = bc->getFixedScalarValue();
            } else if (bc->valueType == BCValueType::VECTOR) {
                // Determine component from scalar field name
                if (phi.name == "Ux" || phi.name == "U_x")
                {
                    boundaryValue = bc->vectorValue.x;
                }
                else if (phi.name == "Uy" || phi.name == "U_y")
                {
                    boundaryValue = bc->vectorValue.y;
                }
                else if (phi.name == "Uz" || phi.name == "U_z") 
                {
                    boundaryValue = bc->vectorValue.z;
                }
                else 
                {
                    // Fallback: assume zero normal gradient
                    return cellGradient;
                }
            }
            else
            {
                // Unknown BC value type: fallback to cell gradient
                return cellGradient;
            }

            Scalar cellValue = phi[face.ownerCell];
            
            // Calculate normal distance from cell center to face
            Scalar d_n = dot(face.d_Pf, face.normal);

            if (std::abs(d_n) < GRADIENT_TOLERANCE)
            {
                return cellGradient;
            }
            
            // Calculate normal gradient: ∂φ/∂n = (φ_boundary - φ_cell) / d_n
            Scalar normalGradient = (boundaryValue - cellValue) / d_n;
            
            // Project cell gradient onto tangential directions
            Vector tangentialGradient = 
                cellGradient - dot(cellGradient, face.normal) * face.normal;
            
            return tangentialGradient + normalGradient * face.normal;
        }
        
        case BCType::ZERO_GRADIENT:
        {
            return cellGradient;
        }
            
        case BCType::FIXED_GRADIENT:
        {
            // For fixed gradient BC, apply the specified gradient value
            Scalar specifiedGradient = bc->getFixedScalarGradient();
            
            // Project cell gradient onto tangential directions and combine
            // with specified normal gradient
            Vector tangentialGradient = 
                cellGradient - dot(cellGradient, face.normal) * face.normal;
            
            return tangentialGradient + specifiedGradient * face.normal;
        }
        
        case BCType::SYMMETRY:
        {
            // Symmetry: zero normal gradient for scalars
            // Remove normal component of gradient
            return cellGradient - dot(cellGradient, face.normal) * face.normal;
        }

        default:
            // For unknown/unhandled BC types, assume zero-gradient
            return cellGradient;
    }
}