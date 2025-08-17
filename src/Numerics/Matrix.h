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
 * @brief Time discretization schemes for transient simulations
 */
enum class TimeScheme
{
    Steady,     ///< Steady-state problems
    Transient   ///< Transient problems
};

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
    /// References to mesh data and numerical schemes
    const std::vector<Face>& allFaces;
    const std::vector<Cell>& allCells;
    const BoundaryConditions& bcManager;
    const GradientScheme& gradScheme;

    /// Cell-centered gradient fields
    VectorField gradP;
    VectorField gradUx;
    VectorField gradUy;
    VectorField gradUz;
    VectorField gradk;
    VectorField gradOmega;

    /// Face-level interpolated data
    FaceFluxField mdotFaces;
    FaceVectorField gradP_f;
    FaceVectorField gradUx_f;
    FaceVectorField gradUy_f;
    FaceVectorField gradUz_f;
    FaceVectorField gradk_f;
    FaceVectorField gradOmega_f;

    /// Cache validity flag for current iteration
    bool cachesValid = false;

    /// Linear system matrix A and right-hand side vector b
    Eigen::SparseMatrix<Scalar> A_matrix;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> b_vector;

    /**
     * @brief Constructor for matrix assembly
     * @param faces Reference to face data
     * @param cells Reference to cell data
     * @param boundaryConds Reference to boundary conditions
     * @param gradientScheme Reference to gradient scheme
     */
    Matrix
    (
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells,
        const BoundaryConditions& boundaryConds,
        const GradientScheme& gradientScheme
    );

    /**
     * @brief Refresh iteration caches for current solution
     * @param p Current pressure field
     * @param U Current velocity field
     * @param rho Fluid density
     * @param turbulenceModel Pointer to turbulence model (optional)
     */
    void refreshIterationCaches
    (
        const PressureField& p,
        const VelocityField& U,
        Scalar rho,
        const KOmegaSST* turbulenceModel
    );

    /**
     * @brief Build momentum equation matrix
     * @param fieldName Name of the velocity component field
     * @param phi Current field values
     * @param phi_old Previous time step values
     * @param phi_source Source term values
     * @param rho Fluid density
     * @param Gamma Diffusion coefficient field
     * @param timeScheme Time discretization scheme
     * @param dt Time step size
     * @param theta Time integration parameter
     * @param grad_phi Cell-centered gradients
     * @param grad_phi_f Face-centered gradients
     * @param convScheme Convection scheme
     */
    void buildMomentumMatrix
    (
        const std::string& fieldName,
        const ScalarField& phi,
        const ScalarField& phi_old,
        const ScalarField& phi_source,
        Scalar rho,
        const ScalarField& Gamma,
        TimeScheme timeScheme,
        Scalar dt,
        Scalar theta,
        const VectorField& grad_phi,
        const FaceVectorField& grad_phi_f,
        const ConvectionScheme& convScheme
    );

    /**
     * @brief Build pressure correction matrix
     * @param massFlux Face mass fluxes
     * @param a_Ux Momentum matrix diagonal for x-velocity
     * @param a_Uy Momentum matrix diagonal for y-velocity
     * @param a_Uz Momentum matrix diagonal for z-velocity
     * @param rho Fluid density
     */
    void buildPressureMatrix
    (
        const FaceFluxField& massFlux,
        const ScalarField& a_Ux,
        const ScalarField& a_Uy,
        const ScalarField& a_Uz,
        Scalar rho
    );

    /**
     * @brief Build scalar transport matrix
     * @param fieldName Name of the scalar field
     * @param phi Current field values
     * @param phi_old Previous time step values
     * @param U_field Velocity field for convection
     * @param phi_source Source term values
     * @param rho Fluid density
     * @param Gamma Diffusion coefficient field
     * @param timeScheme Time discretization scheme
     * @param dt Time step size
     * @param theta Time integration parameter
     * @param convScheme Convection scheme
     */
    void buildScalarTransportMatrix
    (
        const std::string& fieldName,
        const ScalarField& phi,
        const ScalarField& phi_old,
        const VectorField& U_field,
        const ScalarField& phi_source,
        Scalar rho,
        const ScalarField& Gamma,
        TimeScheme timeScheme,
        Scalar dt,
        Scalar theta,
        const ConvectionScheme& convScheme
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
     * @brief Set face mass fluxes for current iteration
     * @param mdot_in Input face mass flux field
     */
    void setFaceMassFluxes(const FaceFluxField& mdot_in) 
    {
        mdotFaces = mdot_in; 
    }

private:
    /// Storage for matrix assembly
    std::vector<Eigen::Triplet<Scalar>> tripletList;
    /// Mapping from face indices to boundary patches
    std::map<size_t, const BoundaryPatch*> faceToPatchMap;

    /**
     * @brief Clear matrix and vector storage
     */
    void clear();

    /**
     * @brief Resolve field name for boundary condition lookup
     * @param fieldName Input field name
     * @return Resolved field name for BC matching
     */
    std::string resolveBCFieldName(const std::string& fieldName) const;

    /**
     * @brief Compute Dirichlet boundary condition value
     * @param bc Boundary condition data
     * @param fieldName Field name for component extraction
     * @return Computed boundary value
     */
    Scalar computeDirichletValue
    (
        const BoundaryData* bc,
        const std::string& fieldName
    ) const;
};

#endif