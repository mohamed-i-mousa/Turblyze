#ifndef MATRIXCONSTRUCTOR_H
#define MATRIXCONSTRUCTOR_H

#include <vector>
#include <string>
#include <map>

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>

#include "Cell.h"
#include "Face.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "ConvectionScheme.h"
#include "GradientScheme.h"

class KOmegaSST;

/**
 * @brief Matrix assembly class for finite volume discretization
 * 
 * This class handles the assembly of linear systems Ax=b for various
 * transport equations including momentum, pressure correction, and
 * scalar transport. It manages gradient fields, face-level data,
 * and boundary condition application.
 */
class Matrix
{
public:
    /**
     * @brief Constructor for matrix assembly
     * @param faces Reference to face data
     * @param cells Reference to cell data
     * @param boundaryConds Reference to boundary conditions
     */
    Matrix
    (
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells,
        const BoundaryConditions& boundaryConds
    );

    /**
     * @brief Deprecated - no longer using caches
     */
    void refreshIterationCaches
    (
        const PressureField& p,
        const VelocityField& U,
        const KOmegaSST* turbulenceModel
    );

    /**
     * @brief Build momentum equation matrix
     * @param phi Current field values
     * @param phi_source Source term values
     * @param Gamma Diffusion coefficient field
     * @param convScheme Convection scheme
     * @param gradScheme Gradient scheme
     * @param fieldName Field name
     */
    void buildMatrix
    (
        const ScalarField& phi,
        const ScalarField& phi_source,
        const VectorField& U_field,
        const ScalarField& Gamma,
        const ConvectionScheme& convScheme,
        const GradientScheme& gradScheme,
        const std::string& fieldName
    );

    /**
     * @brief Build pressure correction matrix
     * @param RhieChowFlowRate Face volume flow rate 
     * @param D_f Face-based diffusion coefficients
     */
    void buildPressureCorrectionMatrix
    (
        const FaceFluxField& RhieChowFlowRate,
        const FaceFluxField& D_f
    );

    /**
     * @brief Get the assembled matrix A
     * @return Reference to the sparse matrix A
     */
    const Eigen::SparseMatrix<Scalar>& getMatrixA() const 
    {
        return A_matrix; 
    }

    /**
     * @brief Get the right-hand side vector b
     * @return Reference to the right-hand side vector
     */
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& getVectorB() const
    {
        return b_vector;
    }

    /**
     * @brief Apply implicit under-relaxation to the linear system
     * @param alpha Relaxation factor (0 < alpha < 1)
     * @param phi_prev Previous iteration values
     */
    void relax(Scalar alpha, const ScalarField& phi_prev);

    /**
     * @brief Set face fluxes for current iteration
     * @param mDot Input face flow rate field
     */

private:
    /// References to mesh data and numerical schemes
    const std::vector<Face>& allFaces;
    const std::vector<Cell>& allCells;
    const BoundaryConditions& bcManager;


    /// Linear system matrix A and right-hand side vector b
    Eigen::SparseMatrix<Scalar> A_matrix;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> b_vector;
    /// Storage for matrix assembly
    std::vector<Eigen::Triplet<Scalar>> tripletList;
    /// Mapping from face indices to boundary patches
    std::map<size_t, const BoundaryPatch*> faceToPatchMap;

    /**
     * @brief Clear matrix and vector storage
     */
    void clear();
    
    /**
     * @brief Get gradient field for specified field type
     * @param fieldName Field name to identify appropriate gradient
     * @param phi Scalar field (for non-cached fields)
     * @param[out] grad_phi_f Face gradient field
     * @return Pointer to appropriate cell-centered gradient field
     */
    const VectorField* getGradientField
    (
        const std::string& fieldName,
        const ScalarField& phi,
        FaceVectorField& grad_phi_f
    );
};

#endif