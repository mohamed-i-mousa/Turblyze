/******************************************************************************
 * @file GradientScheme.hpp
 * @brief Gradient computation schemes for finite volume discretization
 *
 * @details This header defines the GradientScheme class, which provides 
 * robust gradient reconstruction methods for scalar and vector fields on
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
     * @param fieldName Name of the field for BC lookup
     * @param phi Scalar field for gradient calculation
     * @param boundaryFaceValues Optional pre-computed boundary face
     *        values that override BC lookup when provided
     * @return Gradient vector at the specified cell
     */
    Vector cellGradient
    (
        const std::string& fieldName,
        const ScalarField& phi,
        size_t cellIndex,
        const FaceData<Scalar>* boundaryFaceValues = nullptr
    ) const;

    /**
     * @brief Interpolate gradient at a single face
     * 
     * @details
     * For internal faces, computes a corrected face gradient using
     * distance-weighted interpolation of cell gradients with an
     * orthogonal correction to ensure consistency with the direct
     * cell-to-cell difference.
     * 
     * @param fieldName Name of the field for BC lookup
     * @param phi Cell-centered scalar field
     * @param gradPhi_P Gradient at the owner cell
     * @param gradPhi_N Gradient at the neighbor cell
     * @param faceIndex Index of the face
     * @return Gradient vector at the specified face
     */
    Vector faceGradient
    (
        const std::string& fieldName,
        const ScalarField& phi,
        const Vector& gradPhi_P,
        const Vector& gradPhi_N,
        const size_t faceIndex
    ) const;

    /**
     * @brief Apply Barth-Jespersen cell-based gradient limiter
     *
     * @details Scales each cell's gradient so that face-extrapolated
     * values do not exceed the range of neighboring cell values.
     * Prevents overshoot from steep gradients (e.g. near wall cells).
     *
     * @param phi Cell-centered scalar field
     * @param gradPhi Cell-centered gradient field (modified in place)
     */
    void limitGradient
    (
        const ScalarField& phi,
        VectorField& gradPhi
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
        const Vector& gradPhiP,
        const Vector& gradPhiN
    ) const;

    /**
     * @brief Calculate boundary face gradient based on BC type
     * @param face Boundary face
     * @param fieldName Name of the field for BC lookup
     * @param phi Scalar field values
     * @param cellGradient Gradient at the owner cell
     * @return Gradient vector at the boundary face
     */
    Vector calculateBoundaryFaceGradient
    (
        const std::string& fieldName,
        const ScalarField& phi,
        const Vector& cellGradient,
        const Face& face
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
