/******************************************************************************
 * @file Matrix.hpp
 * @brief Matrix assembly for finite volume discretization
 *
 * This header defines the Matrix class, which handles the assembly of sparse
 * linear systems (Ax=b) for transport equations.
 * It manages gradient computation, face-level interpolation, convection and
 * diffusion discretization, and boundary condition application.
 *
 * @class Matrix
 *
 * The Matrix class provides:
 * - Assembly for momentum, pressure correction, and scalar equations
 * - Eigen sparse matrix for efficient storage and solution
 * - Deferred correction for higher-order convection schemes
 * - Non-orthogonal mesh corrections
 * - Implicit under-relaxation (Patankar) for solution stability
 * - Boundary condition integration during assembly
 *****************************************************************************/

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <string>
#include <map>

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>

#include "Cell.hpp"
#include "Face.hpp"
#include "BoundaryConditions.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"
#include "ConvectionScheme.hpp"
#include "GradientScheme.hpp"


class Matrix
{
public:

    /**
     * @brief Constructor for matrix assembly
     * @param faces Reference to face data
     * @param cells Reference to cell data
     * @param boundaryConds Boundary conditions manager
     */
    Matrix
    (
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells,
        const BoundaryConditions& boundaryConds
    );

    /**
     * @brief Build transport equation matrix with convection-diffusion
     * @param phi Current field values
     * @param phiSource Source term field
     * @param flowRateFace Face volumetric flow rates
     * @param Gamma Diffusion coefficient field
     * @param convScheme Convection discretization scheme
     * @param gradScheme Gradient reconstruction scheme
     * @param gradPhi Pre-computed cell gradients of phi
     * @param fieldName Name of the field being solved
     *        (e.g., "Ux", "Uy", "Uz", "k", "omega")
     */
    void buildMatrix
    (
        const ScalarField& phi,
        const ScalarField& phiSource,
        const FaceFluxField& flowRateFace,
        const ScalarField& Gamma,
        const ConvectionScheme& convScheme,
        const GradientScheme& gradScheme,
        const VectorField& gradPhi,
        const std::string& fieldName,
        const FaceData<Scalar>* boundaryFaceValues = nullptr
    );

    /**
     * @brief Build pressure correction matrix for SIMPLE algorithm
     * @param RhieChowFlowRate Face flow rate from Rhie-Chow
     * @param DUf Face diffusion coefficients (D_f = V / a_P)
     * @param pCorr Pressure correction field
     * @param gradScheme Gradient reconstruction scheme
     * @param gradpCorr Pre-computed cell gradients of pCorr
     */
    void buildPressureCorrectionMatrix
    (
        const FaceFluxField& RhieChowFlowRate,
        const FaceFluxField& DUf,
        const ScalarField& pCorr,
        const GradientScheme& gradScheme,
        const VectorField& gradpCorr
    );

// Accessor methods

    /**
     * @brief Get assembled sparse matrix A (const)
     * @return Const reference to coefficient matrix
     */
    const Eigen::SparseMatrix<Scalar>& getMatrixA() const
    {
        return matrixA_;
    }

    /**
     * @brief Get assembled sparse matrix A (non-const)
     * @return Mutable reference to coefficient matrix
     */
    Eigen::SparseMatrix<Scalar>& getMatrixA()
    {
        return matrixA_;
    }

    /**
     * @brief Get right-hand side vector b (const)
     * @return Const reference to RHS vector
     */
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&
    getVectorB() const
    {
        return vectorB_;
    }

    /**
     * @brief Get right-hand side vector b (non-const)
     * @return Mutable reference to RHS vector
     */
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& getVectorB()
    {
        return vectorB_;
    }

    /**
     * @brief Apply Patankar implicit under-relaxation
     * @param alpha Relaxation factor (0 < alpha <= 1)
     * @param phiPrev Previous iteration field values
     *
     * Modifies the assembled system (Ax=b):
     * - Diagonal: a_P <- a_P / alpha
     * - RHS: b <- b + ((1-alpha)/alpha) * aPOriginal * phiPrev
     */
    void relax(Scalar alpha, const ScalarField& phiPrev);

private:

// Private members

    /// Mesh and boundary data references
    const std::vector<Face>& allFaces_;
    const std::vector<Cell>& allCells_;
    const BoundaryConditions& bcManager_;

    /// Sparse linear system components
    Eigen::SparseMatrix<Scalar> matrixA_;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vectorB_;

    /// Triplet storage for efficient sparse matrix assembly
    std::vector<Eigen::Triplet<Scalar>> tripletList_;

    /// Cached face counts for triplet list reservation
    size_t numInternalFaces_;
    size_t numBoundaryFaces_;

// Private methods

    /**
     * @brief Clear matrix and vector for new assembly
     */
    void clear();

    /**
     * @brief Reserve triplet list capacity based on mesh topology
     */
    void reserveTripletList();

    /**
     * @brief Assemble internal face for transport equation
     * @param face Internal face to process
     * @param phi Current field values
     * @param flowRateFace Face volumetric flow rates
     * @param Gamma Diffusion coefficient field
     * @param convScheme Convection discretization scheme
     * @param gradScheme Gradient reconstruction scheme
     * @param gradPhi Pre-computed cell gradients
     * @param fieldName Name of the field being solved
     */
    void assembleInternalFace
    (
        const Face& face,
        const ScalarField& phi,
        const FaceFluxField& flowRateFace,
        const ScalarField& Gamma,
        const ConvectionScheme& convScheme,
        const GradientScheme& gradScheme,
        const VectorField& gradPhi,
        const std::string& fieldName
    );

    /**
     * @brief Assemble boundary face for transport equation
     * @param face Boundary face to process
     * @param phi Current field values
     * @param flowRateFace Face volumetric flow rates
     * @param Gamma Diffusion coefficient field
     * @param gradScheme Gradient reconstruction scheme
     * @param gradPhi Pre-computed cell gradients
     * @param fieldName Name of the field being solved
     */
    void assembleBoundaryFace
    (
        const Face& face,
        const ScalarField& phi,
        const FaceFluxField& flowRateFace,
        const ScalarField& Gamma,
        const GradientScheme& gradScheme,
        const VectorField& gradPhi,
        const std::string& fieldName,
        const FaceData<Scalar>* boundaryFaceValues = nullptr
    );

    /**
     * @brief Assemble internal face for pressure correction
     * @param face Internal face to process
     * @param DUf Face diffusion coefficients
     * @param pCorr Pressure correction field
     * @param gradScheme Gradient reconstruction scheme
     * @param gradpCorr Pre-computed pCorr gradients
     */
    void assemblePCorrInternalFace
    (
        const Face& face,
        const FaceFluxField& DUf,
        const ScalarField& pCorr,
        const GradientScheme& gradScheme,
        const VectorField& gradpCorr
    );

    /**
     * @brief Assemble boundary face for pressure correction
     * @param face Boundary face to process
     * @param DUf Face diffusion coefficients
     * @param grad_pCorr Pre-computed pCorr gradients
     */
    void assemblePCorrBoundaryFace
    (
        const Face& face,
        const FaceFluxField& DUf,
        const VectorField& gradpCorr
    );

    /**
     * @brief Extract scalar boundary value from BC data
     * @param bc Boundary data for the patch/field
     * @param fieldName Field name for vector component dispatch
     * @return Scalar boundary value
     */
    static Scalar extractBoundaryScalar
    (
        const BoundaryData& bc,
        const std::string& fieldName
    );
};

#endif // MATRIX_HPP
