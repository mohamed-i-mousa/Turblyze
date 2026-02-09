/******************************************************************************
 * @file GradientScheme.hpp
 * @brief Gradient computation schemes for finite volume discretization
 *
 * This header defines the GradientScheme class, which provides robust
 * gradient reconstruction methods for scalar and vector fields on
 * unstructured finite volume meshes. The implementation uses a weighted
 * least-squares approach with inverse distance squared weighting.
 *
 * @class GradientScheme
 *
 * The GradientScheme class provides:
 * - Cell-centered gradient computation via weighted least-squares
 * - Face-centered gradient interpolation with orthogonal correction
 * - Boundary condition handling for gradient stencils
 * - Distance-weighted averaging for internal face gradients
 *****************************************************************************/

#ifndef GRADIENT_SCHEME_HPP
#define GRADIENT_SCHEME_HPP

#include <vector>
#include <string>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "BoundaryConditions.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"

class GradientScheme
{
public:

    /**
     * @brief Construct gradient scheme with mesh context
     * @param faces Reference to all mesh faces
     * @param cells Reference to all mesh cells
     * @param bc Reference to boundary conditions manager
     */
    GradientScheme
    (
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells,
        const BoundaryConditions& bc
    );

    /**
     * @brief Calculate gradient at a single cell using least-squares.
     * @param cellIndex Index of the cell to compute gradient for
     * @param phi Scalar field for gradient calculation
     * @param fieldName Name of the field for BC lookup
     * @return Gradient vector at the specified cell
     * 
     * Uses weighted least-squares reconstruction including both
     * internal neighbor cells and boundary faces.
     */
    Vector cellGradient
    (
        size_t cellIndex,
        const ScalarField& phi,
        const std::string& fieldName
    ) const;

    /**
     * @brief Interpolate gradient at a single face
     *
     * For internal faces, computes a corrected face gradient using
     * distance-weighted interpolation of cell gradients with an
     * orthogonal correction to ensure consistency with the direct
     * cell-to-cell difference.
     *
     * @param faceIndex Index of the face
     * @param gradPhi_P Gradient at the owner cell
     * @param gradPhi_N Gradient at the neighbor cell
     * @param phi Cell-centered scalar field
     * @param fieldName Name of the field for BC lookup
     * @return Gradient vector at the specified face
     */
    Vector faceGradient
    (
        const size_t faceIndex,
        const Vector& gradPhi_P,
        const Vector& gradPhi_N,
        const ScalarField& phi,
        const std::string& fieldName
    ) const;


private:

// Private methods

    /**
     * @brief Distance-weighted linear interpolation of gradients
     * @param face Internal face for interpolation
     * @param gradPhi_P Gradient at the owner cell
     * @param gradPhi_N Gradient at the neighbor cell
     * @return Linearly interpolated gradient vector at the face
     */
    Vector averageFaceGradient
    (
        const Face& face,
        const Vector& gradPhi_P,
        const Vector& gradPhi_N
    ) const;

    /**
     * @brief Calculate boundary face gradient based on BC type
     * @param face Boundary face
     * @param cellGradient Gradient at the owner cell
     * @param phi Scalar field values
     * @param fieldName Name of the field for BC lookup
     * @return Gradient vector at the boundary face
     */
    Vector calculateBoundaryFaceGradient
    (
        const Face& face,
        const Vector& cellGradient,
        const ScalarField& phi,
        const std::string& fieldName
    ) const;

// Private members

    /// Reference to all mesh faces
    const std::vector<Face>& allFaces_;

    /// Reference to all mesh cells
    const std::vector<Cell>& allCells_;
    
    /// Reference to boundary conditions manager
    const BoundaryConditions& bcManager_;
};

#endif // GRADIENT_SCHEME_HPP