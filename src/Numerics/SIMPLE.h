#ifndef SIMPLE_H
#define SIMPLE_H

#include <vector>
#include <memory>

#include "Face.h"
#include "Cell.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "FaceData.h"
#include "GradientScheme.h"
#include "ConvectionScheme.h"
#include "Matrix.h"
#include "LinearSolvers.h"
#include "KOmegaSST.h"


/**
 * @brief SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) solver
 * 
 * This class implements the SIMPLE algorithm for solving incompressible
 * Navier-Stokes equations using finite volume discretization. It handles
 * velocity-pressure coupling through pressure correction and includes
 * support for turbulence modeling.
 */
class SIMPLE {
public:
    /**
     * @brief Constructor for SIMPLE solver
     * @param faces Reference to face data
     * @param cells Reference to cell data
     * @param bc Reference to boundary conditions
     * @param gradScheme Reference to gradient scheme
     * @param convScheme Reference to convection scheme
     */
    SIMPLE
    (
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells, 
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionScheme& convScheme
    );

    /**
     * @brief Main SIMPLE algorithm execution
     */
    void solve();
    
    /**
     * @brief Solve momentum equations for velocity components
     */
    void solveMomentumEquations();
    
    /**
     * @brief Calculate mass fluxes using Rhie-Chow interpolation
     */
    void calculateRhieChowMassFlux();

    /**
     * @brief Solve pressure correction equation
     */
    void solvePressureCorrection();
    
    /**
     * @brief Correct velocity field using pressure correction
     */
    void correctVelocity();
    
    /**
     * @brief Correct pressure field using pressure correction
     */
    void correctPressure();
    
    /**
     * @brief Correct mass fluxes using updated velocity field
     */
    void correctMassFluxes();
    
    /**
     * @brief Check if solution has converged
     * @return true if converged, false otherwise
     */
    bool checkConvergence();

    /**
     * @brief Get velocity field
     * @return Reference to velocity field
     */
    const VectorField& getVelocity() const { return U; }
    
    /**
     * @brief Get pressure field
     * @return Reference to pressure field
     */
    const ScalarField& getPressure() const { return p; }
    
    /**
     * @brief Get mass flux field
     * @return Reference to mass flux field
     */
    const FaceFluxField& getRhieChowMassFlux() const 
    { 
        return RhieChowMassFlux; 
    }
    
    /**
     * @brief Get x-component of velocity field
     * @return Scalar field of x-velocity components
     */
    ScalarField getVelocityX() const;
    
    /**
     * @brief Get y-component of velocity field
     * @return Scalar field of y-velocity components
     */
    ScalarField getVelocityY() const;
    
    /**
     * @brief Get z-component of velocity field
     * @return Scalar field of z-velocity components
     */
    ScalarField getVelocityZ() const;

    /**
     * @brief Set under-relaxation factors
     * @param alpha_U Velocity under-relaxation factor
     * @param alpha_p Pressure under-relaxation factor
     */
    void setRelaxationFactors(Scalar alpha_U, Scalar alpha_p);
    
    /**
     * @brief Set convergence tolerance
     * @param tol Convergence tolerance value
     */
    void setConvergenceTolerance(Scalar tol);
    
    /**
     * @brief Set maximum number of iterations
     * @param maxIter Maximum iteration count
     */
    void setMaxIterations(int maxIter);
    
    /**
     * @brief Enable or disable turbulence modeling
     * @param enable Whether to enable turbulence modeling
     */
    void enableTurbulenceModeling(bool enable = true);
    
    /**
     * @brief Set physical properties
     * @param rho_new Fluid density
     * @param mu_new Dynamic viscosity
     */
    void setPhysicalProperties(Scalar rho_new, Scalar mu_new) 
    {
        this->rho = rho_new;
        this->mu = mu_new;
        this->nu = mu_new / rho_new;
    }
    
    /**
     * @brief Get turbulent kinetic energy field
     * @return Pointer to k field (null if turbulence disabled)
     */
    const ScalarField* getTurbulentKineticEnergy() const;
    
    /**
     * @brief Get specific dissipation rate field
     * @return Pointer to omega field (null if turbulence disabled)
     */
    const ScalarField* getSpecificDissipationRate() const;
    
    /**
     * @brief Get turbulent viscosity field
     * @return Pointer to nu_t field (null if turbulence disabled)
     */
    const ScalarField* getTurbulentViscosity() const;
    
    /**
     * @brief Get wall distance field
     * @return Pointer to wall distance field (null if turbulence disabled)
     */
    const ScalarField* getWallDistance() const;

private:
    /// Mesh references
    const std::vector<Face>& allFaces;
    const std::vector<Cell>& allCells;
    const BoundaryConditions& bcManager;
    const GradientScheme& gradientScheme;
    const ConvectionScheme& convectionScheme;

    /// Physical properties
    Scalar rho;           ///< Fluid density
    Scalar mu;            ///< Dynamic viscosity
    Scalar nu;            ///< Kinematic viscosity

    /// Algorithm parameters
    Scalar alpha_U;         ///< Under-relaxation factor for velocity
    Scalar alpha_p;         ///< Under-relaxation factor for pressure
    int maxIterations;      ///< Maximum number of iterations
    Scalar tolerance;       ///< Convergence tolerance
    bool enableTurbulence;  ///< Enable turbulence modeling
    
    /// Turbulence model
    std::unique_ptr<KOmegaSST> turbulenceModel;

    /// Solution fields
    VectorField U;        ///< Velocity field
    ScalarField p;        ///< Pressure field
    ScalarField p_prime;  ///< Pressure correction field
    Scalar lastPressureCorrectionRMS = S(1e9); ///< Track p' RMS before reset
    
    /// Face-based fields for Rhie-Chow interpolation
    FaceFluxField RhieChowMassFlux;      ///< Mass flux through faces
    FaceVectorField U_face;              ///< Face velocity field

    /// Previous-iteration fields (for under-relaxation effects at faces)
    VectorField U_prev;          ///< Cell-centered velocity from previous iteration
    FaceVectorField U_face_prev; ///< Face velocity from previous iteration 

    /// Momentum equation coefficients per component (needed for Rhie-Chow)
    ScalarField a_Ux;     ///< Diagonal coefficients of U_x momentum equation
    ScalarField a_Uy;     ///< Diagonal coefficients of U_y momentum equation
    ScalarField a_Uz;     ///< Diagonal coefficients of U_z momentum equation
    VectorField H_U;      ///< H/A terms for momentum equations
    
    /// Gradient fields
    VectorField gradP;    ///< Pressure gradient
    
    /// Matrix constructor and solver objects
    std::unique_ptr<Matrix> matrixConstruct;
    
    /**
     * @brief Initialize solution fields
     * @param initialVelocity Initial velocity field
     * @param initialPressure Initial pressure field
     */
    void initialize(const Vector& initialVelocity, Scalar initialPressure);
    
    /**
     * @brief Calculate mass imbalance across domain
     * @return Mass imbalance value
     */
    Scalar calculateMassImbalance() const;
    
    /**
     * @brief Calculate velocity residual
     * @return Velocity residual value
     */
    Scalar calculateVelocityResidual() const;
    
    /**
     * @brief Calculate pressure residual
     * @return Pressure residual value
     */
    Scalar calculatePressureResidual() const;
};

#endif