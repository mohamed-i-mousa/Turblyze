/******************************************************************************
 * @file Matrix.hpp
 * @brief Matrix assembly for finite volume discretization
 *
 * @details This header defines: TransportEquation struct and the Matrix class.
 * TransportEquation bundles all data describing a scalar transport equation
 * (field, convection, diffusion, source, gradients). The Matrix class
 * handles assembly of sparse linear systems (Ax=b) for any transport
 * equation type — momentum, pressure correction, and turbulence.
 *
 * @class Matrix
 *
 * The Matrix class provides:
 * - Unified assembly for all transport equation types
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
#include <optional>
#include <functional>

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>

#include "Cell.hpp"
#include "Face.hpp"
#include "BoundaryConditions.hpp"
#include "CellData.hpp"
#include "FaceData.hpp"
#include "ConvectionScheme.hpp"
#include "GradientScheme.hpp"


/// Non-owning optional reference (absent = std::nullopt)
template<typename T>
using OptionalRef = std::optional<std::reference_wrapper<const T>>;


/**
 * @brief Data describing a scalar transport equation
 *
 * @details 
 * - Combines the field value, convection, diffusion, source,
 *   and gradient data needed by Matrix::buildMatrix().
 * 
 * - Organised by physics terms:
 *     field, transient (placeholder), convection, diffusion,
 *     source, gradient reconstruction, boundary overrides.
 */
struct TransportEquation
{

// Field

    /// Name of the field ("Ux", "k", "pCorr", etc.)
    std::string fieldName;

    /// Current cell-centered field values (mutable for zero-copy solve)
    ScalarField& phi;

// Transient (placeholder for future)

// Convection: div(F * phi)

    /// Face volumetric flow rates (nullopt = no convection)
    OptionalRef<FaceFluxField> flowRate = std::nullopt;

    /// Convection discretization scheme (nullopt = no convection)
    OptionalRef<ConvectionScheme> convScheme = std::nullopt;

// Diffusion: div(Gamma * grad(phi))

    /// Cell-centered diffusion coefficient
    OptionalRef<ScalarField> Gamma = std::nullopt;

    /// Pre-interpolated face diffusion coefficient
    OptionalRef<FaceFluxField> GammaFace = std::nullopt;

// Source

    /// Explicit source term field
    const ScalarField& source;

// Gradient reconstruction

    /// Pre-computed cell gradients of phi
    const VectorField& gradPhi;

    /// Gradient reconstruction scheme
    const GradientScheme& gradScheme;

// Boundary overrides

    /// Dynamic boundary face values (e.g. omega wall-function values)
    OptionalRef<FaceData<Scalar>> boundaryFaceValues = std::nullopt;
};


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
     * @brief Build transport equation matrix
     * @param equation Transport equation data
     */
    void buildMatrix(const TransportEquation& equation);

// Accessor methods

    /**
     * @brief Get assembled sparse matrix A (const)
     * @return Const reference to coefficient matrix
     */
    const Eigen::SparseMatrix<Scalar>& matrixA() const
    {
        return matrixA_;
    }

    /**
     * @brief Get assembled sparse matrix A (non-const)
     * @return Mutable reference to coefficient matrix
     */
    Eigen::SparseMatrix<Scalar>& matrixA()
    {
        return matrixA_;
    }

    /**
     * @brief Get right-hand side vector b (const)
     * @return Const reference to RHS vector
     */
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&
    vectorB() const
    {
        return vectorB_;
    }

    /**
     * @brief Get right-hand side vector b (non-const)
     * @return Mutable reference to RHS vector
     */
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& 
    vectorB()
    {
        return vectorB_;
    }

    /**
     * @brief Apply Patankar implicit under-relaxation
     * @param alpha Relaxation factor (0 < alpha <= 1)
     * @param phiPrev Previous iteration field values
     *
     * @details 
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
     * @brief Assemble internal face contributions
     * @param face Internal face to process
     * @param equation Transport equation data
     */
    void assembleInternalFace
    (
        const Face& face,
        const TransportEquation& equation
    );

    /**
     * @brief Assemble boundary face contributions
     * @param face Boundary face to process
     * @param equation Transport equation data
     */
    void assembleBoundaryFace
    (
        const Face& face,
        const TransportEquation& equation
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
