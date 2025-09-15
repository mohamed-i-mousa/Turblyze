/******************************************************************************
 * @file GradientScheme.h
 * @brief Gradient computation schemes for finite volume discretization
 * 
 * This class provides robust gradient reconstruction methods for scalar and
 * vector fields on unstructured finite volume meshes. The implementation
 * uses weighted least-squares approach with inverse distance squared 
 * weighting and includes gradient limiting for stability. Both cell-centered
 * and face-centered gradients are computed with proper boundary condition 
 * handling.
 * 
 * @class GradientScheme
 * 
 * The GradientScheme class offers:
 * - Weighted least-squares gradient reconstruction
 * - Automatic regularization for singular matrices  
 * - Barth-Jesperson gradient limiting for boundedness
 * - Face gradient interpolation with orthogonal corrections
 * - Seamless boundary condition integration
 * - Support for both scalar and vector fields
 * 
 * Key numerical features:
 * - Robust matrix solver with LLT/LU fallback mechanisms
 * - Distance-weighted interpolation for face gradients
 * - Normal/tangential gradient decomposition at boundaries
 * - Consistent treatment of non-orthogonal mesh corrections
 *****************************************************************************/

#ifndef GRADIENTSCHEME_H
#define GRADIENTSCHEME_H

#include <vector>
#include <string>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "BoundaryConditions.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"
#include "linearInterpolation.hpp"

/**
 * @file GradientScheme.h
 * @brief Gradient reconstruction and interpolation schemes
 * 
 * For each cell P, we solve the overdetermined system:
 *   phi_N - phi_P = ∇φ_P · (r_N - r_P) + O(|r_N - r_P|²)
 * 
 * This leads to the normal equations:  A^T A ∇φ_P = A^T b
 * where:
 *   A_ij = w_i * (r_i - r_P)_j
 *   b_i = w_i * (phi_i - phi_P)
 *   w_i = 1/|r_i - r_P|² (inverse distance squared weighting)
 */

/**
 * @brief Element-level gradient reconstruction and interpolation
 * 
 * Provides cell-level gradient computation and face-level gradient
 * interpolation using least-squares reconstruction.
 * 
 */
class GradientScheme
{
public:
    GradientScheme() = default;

    /**
     * @brief Calculate gradient at a single cell using least-squares
     * @param cellIndex Index of the cell to compute gradient for
     * @param phi Scalar field for gradient calculation
     * @param allCells Vector of all cells in the mesh
     * @return Gradient vector at the specified cell
     */
    Vector CellGradient
    (
        size_t cellIndex,
        const ScalarField& phi,
        const std::vector<Cell>& allCells
    ) const;

    /**
     * @brief Interpolate gradient at a single face
     * @param faceIndex Index of the face to interpolate gradient for
     * @param grad_phi Cell-centered gradient field
     * @param phi Cell-centered scalar field
     * @param allCells Vector of all cells in the mesh
     * @param allFaces Vector of all faces in the mesh
     * @param boundaryConditions Boundary conditions manager
     * @param fieldName Name of the field for boundary condition lookup
     * @return Gradient vector at the specified face
     */
    Vector FaceGradient
    (
        const size_t faceIndex,
        const Vector& grad_phi_P,
        const Vector& grad_phi_N,
        const ScalarField& phi,
        const std::vector<Cell>& allCells,
        const std::vector<Face>& allFaces,
        const BoundaryConditions& boundaryConditions,
        const std::string& fieldName
    ) const;

    /**
     * @brief Linear interpolation of gradients at a face
     * @param face Internal face for interpolation
     * @param grad_phi Cell-centered gradient field
     * @return Linearly interpolated gradient vector at the face
     */
    Vector averageFaceGradient
    (
        const Face& face,
        const Vector& grad_phi_P,
        const Vector& grad_phi_N
    ) const;

    
private:
    /**
     * @brief Calculate boundary face gradient based on boundary condition type
     * @param face Boundary face
     * @param cellGradient Gradient at the owner cell
     * @param phi Scalar field values
     * @param boundaryConditions Boundary conditions manager
     * @param fieldName Name of the field for boundary condition lookup
     * @return Gradient vector at the boundary face
     */
    Vector calculateBoundaryFaceGradient
    (
        const Face& face,
        const Vector& cellGradient,
        const ScalarField& phi,
        const BoundaryConditions& boundaryConditions,
        const std::string& fieldName
    ) const;
};


#endif