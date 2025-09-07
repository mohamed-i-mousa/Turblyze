#include "GradientScheme.h"

#include <stdexcept>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <eigen3/Eigen/Dense>


Vector GradientScheme::CellGradient
(
    size_t cellIndex,
    const ScalarField& phi,
    const std::vector<Cell>& allCells
) const
{
    const Cell& cell = allCells[cellIndex];

    Eigen::Matrix<Scalar,3,3> ATA;
    Eigen::Matrix<Scalar,3,1> ATb;
    Eigen::Matrix<Scalar,3,1> r_vector;
    
    ATA.setZero();
    ATb.setZero();

    Scalar totalWeight = S(0.0);
    Scalar maxDistance = 0.0;

    for (size_t neighborId : cell.neighborCellIndices()) 
    {
        if (neighborId >= allCells.size()) continue; 
        
        const Cell& neighbor = allCells[neighborId];
        Vector r = neighbor.centroid() - cell.centroid();
        
        Scalar r_mag_sq = r.magnitudeSquared();
        
        maxDistance = std::max(maxDistance, r.magnitude());

        Scalar w = S(1.0) / (r_mag_sq + vSmallValue);
        totalWeight += w;

        r_vector << r.x(), r.y(), r.z();
        ATA.noalias() += w * (r_vector * r_vector.transpose());
        
        Scalar delta_phi = phi[neighborId] - phi[cellIndex];
        ATb.noalias() += w * delta_phi * r_vector;
    }

    Scalar regularization = totalWeight * 1e-12;
    ATA(0,0) += regularization;
    ATA(1,1) += regularization;
    ATA(2,2) += regularization;

    Eigen::LLT<Eigen::Matrix<Scalar,3,3>> llt(ATA);

    if (llt.info() == Eigen::Success) 
    {
        Eigen::Matrix<Scalar,3,1> g = llt.solve(ATb);
        
        // Check solution validity
        Scalar gradMag = g.norm();
        Scalar phiRange = 0.0;
        for (size_t neighborId : cell.neighborCellIndices()) {
            if (neighborId < allCells.size()) {
                phiRange = std::max(phiRange, 
                    std::abs(phi[neighborId] - phi[cellIndex]));
            }
        }
        
        // Barth-Jesperson type limiter
        if (gradMag * maxDistance > 10.0 * phiRange && phiRange > 1e-10)
        {
            g *= (10.0 * phiRange) / (gradMag * maxDistance);
        }
        
        return Vector(g(0), g(1), g(2));    
    }
    else
    {
        Eigen::FullPivLU<Eigen::Matrix<Scalar,3,3>> lu(ATA);

        if (lu.isInvertible()) 
        {
            Eigen::Matrix<Scalar,3,1> g = lu.solve(ATb);
            return Vector(g(0), g(1), g(2));
        }
        else
        {
            throw std::runtime_error
            (
                "Gradient computation failed!"
            );        
        }
    }
}

Vector GradientScheme::FaceGradient
(
    const size_t faceIndex,
    const Vector& grad_phi_P,
    const Vector& grad_phi_N,
    const ScalarField& phi,
    const std::vector<Cell>& allCells,
    const std::vector<Face>& allFaces,
    const BoundaryConditions& boundaryConditions,
    const std::string& fieldName
) const
{
    const Face& face = allFaces[faceIndex];
    size_t P = face.ownerCell();
    
    if (face.isBoundary())
    {
        return 
            calculateBoundaryFaceGradient
            (
                face, grad_phi_P, phi, boundaryConditions, fieldName
            );
    }
    else
    {
        size_t N = face.neighborCell().value();
        
        Vector d_PN = allCells[N].centroid() - allCells[P].centroid();

        Scalar d_PN_mag = d_PN.magnitude();
        
        Vector e_PN = d_PN / (d_PN_mag + vSmallValue);

        Vector grad_avg = averageFaceGradient(face, grad_phi_P, grad_phi_N);
        
        Scalar phi_diff = phi[N] - phi[P];

        Scalar correction = (phi_diff / d_PN_mag) - dot(grad_avg, e_PN);
        
        return grad_avg + correction * e_PN;
    }
}

Vector GradientScheme::averageFaceGradient
(
    const Face& face,
    const Vector& grad_phi_P,
    const Vector& grad_phi_N
) const
{
    if (face.isBoundary())
    {
        throw std::invalid_argument
        (
            "VectorLinearInterpolation: Cannot interpolate on boundary face"
        );
    }

    Scalar d_Pf = face.d_Pf_mag();
    Scalar d_Nf = face.d_Nf_mag().value();
    Scalar total_dist = d_Pf + d_Nf;

    Scalar g_P = d_Nf / (total_dist + vSmallValue);
    Scalar g_N = d_Pf / (total_dist + vSmallValue);

    return g_P * grad_phi_P + g_N * grad_phi_N;
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
    for (const auto& p : boundaryConditions.patches())
    {
        if (face.id() >= p.firstFaceIndex() && face.id() <= p.lastFaceIndex())
        {
            patch = &p;
            break;
        }
    }

    // If patch not found, use cell gradient
    if (!patch) return cellGradient;

    // Get the boundary condition for this field 
    const BoundaryData* bc = 
        boundaryConditions.fieldBC(patch->patchName(), fieldName);
    
    // If no BC specified, use cell gradient
    if (!bc) return cellGradient;
    
    switch (bc->type()) {
        case BCType::FIXED_VALUE: {
            // For fixed value BC, calculate normal gradient from boundary-cell
            // difference and keep tangential components from cell gradient
            
            Scalar boundaryValue = S(0.0);
            if (bc->valueType() == BCValueType::SCALAR) 
            {
                boundaryValue = bc->fixedScalarValue();
            } 
            else if (bc->valueType() == BCValueType::VECTOR) 
            {
                // Determine component from scalar field name
                if (phi.name == "Ux" || phi.name == "U_x")
                {
                    boundaryValue = bc->vectorValue().x();
                }
                else if (phi.name == "Uy" || phi.name == "U_y")
                {
                    boundaryValue = bc->vectorValue().y();
                }
                else if (phi.name == "Uz" || phi.name == "U_z") 
                {
                    boundaryValue = bc->vectorValue().z();
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

            Scalar cellValue = phi[face.ownerCell()];
            
            // Calculate normal distance from cell center to face
            Scalar d_n = dot(face.d_Pf(), face.normal());
            
            // Calculate normal gradient: ∂φ/∂n = (φ_boundary - φ_cell) / d_n
            Scalar normalGradient = (boundaryValue - cellValue) / (d_n + vSmallValue);
            
            // Project cell gradient onto tangential directions
            Vector tangentialGradient = 
                cellGradient - dot(cellGradient, face.normal()) * face.normal();
            
            return tangentialGradient + normalGradient * face.normal();
        }
        
        case BCType::ZERO_GRADIENT:
        {
            return cellGradient;
        }
            
        case BCType::FIXED_GRADIENT:
        {
            // For fixed gradient BC, apply the specified gradient value
            Scalar specifiedGradient = bc->fixedScalarGradient();
            
            // Project cell gradient onto tangential directions and combine
            // with specified normal gradient
            Vector tangentialGradient = 
                cellGradient - dot(cellGradient, face.normal()) * face.normal();
            
            return tangentialGradient + specifiedGradient * face.normal();
        }

        default:
            // For unknown/unhandled BC types, assume zero-gradient
            return cellGradient;
    }
}