#ifndef KOMEGASST_H
#define KOMEGASST_H

#include <vector>
#include <memory>

#include "Scalar.h"
#include "Vector.h"
#include "Cell.h"
#include "Face.h"
#include "CellData.h"
#include "FaceData.h"
#include "BoundaryConditions.h"
#include "GradientScheme.h"
#include "LinearSolvers.h"

class Matrix;


/* k-omega SST (Shear Stress Transport) Turbulence Model
 * 
 * This class implements the k-omega SST turbulence model developed by Menter.
 * The model combines the robustness of the k-omega model in the near-wall region 
 * with the freestream independence of the k-epsilon model in the far field.
 * 
 * The model solves two transport equations:
 * - Turbulent kinetic energy (k)
 * - Specific dissipation rate (omega)
 * 
 * Key features:
 * - Automatic switching between k-omega and k-epsilon formulations
 * - Enhanced near-wall treatment
 * - Wall distance calculation using Poisson equation
 * - Proper wall shear stress computation
 */
class KOmegaSST {
public:
    // Constructor
    KOmegaSST(const std::vector<Face>& faces,
              const std::vector<Cell>& cells,
              const BoundaryConditions& bc,
              const GradientScheme& gradScheme);
    
    ~KOmegaSST();


    void initialize(const VectorField& U_field, Scalar rho, Scalar mu_lam);

    // New overload that avoids constructing a std::vector every call
    void solve(const VectorField& U_field,
               const VectorField& gradUx,
               const VectorField& gradUy,
               const VectorField& gradUz,
               Scalar rho,
               Scalar mu_lam);

    // Legacy interface kept for internal use
    void solve(const VectorField& U_field,
               const std::vector<VectorField>& gradU,
               Scalar rho,
               Scalar mu_lam);

    void calculateWallDistance();

    void solveOmegaEquation(const VectorField& U_field,
                           const std::vector<VectorField>& gradU,
                           Scalar rho,
                           Scalar mu_lam);


    void applyNearWallTreatmentOmega(Scalar rho, Scalar mu_lam);

    void solveKEquation(const VectorField& U_field,
                       const std::vector<VectorField>& gradU,
                       Scalar rho,
                       Scalar mu_lam);

    void calculateTurbulentViscosity(const VectorField& U_field,
                                    const std::vector<VectorField>& gradU,
                                    Scalar rho,
                                    Scalar mu_lam);

    /**
     * @brief Apply wall corrections for turbulent viscosity
     * Ensures μₜ = 0 at walls and proper behavior in viscous sublayer
     */
    void applyWallCorrections();

    /**
     * @brief Calculate wall shear stress using parallel velocity only
     * τ_wall = μ_eff * (∂U_parallel/∂n)_wall
     * @param U_field Current velocity field
     * @param mu_lam Laminar dynamic viscosity
     */
    void calculateWallShearStress(const VectorField& U_field, Scalar mu_lam);

    // Getters for turbulence fields
    const ScalarField& getK() const { return k; }
    const ScalarField& getOmega() const { return omega; }
    const ScalarField& getTurbulentViscosity() const { return mu_t; }
    const ScalarField& getWallDistance() const { return wallDistance; }
    const ScalarField& getWallShearStress() const { return wallShearStress; }

    /**
     * @brief Get effective viscosity (laminar + turbulent)
     * @param mu_lam Laminar dynamic viscosity
     * @return Effective viscosity field
     */
    ScalarField getEffectiveViscosity(Scalar mu_lam) const;

    // Setters for model parameters
    void setModelConstants(Scalar C_mu = 0.09, Scalar beta_1 = 0.075, Scalar beta_2 = 0.0828);
    
    // Transient time stepping parameters
    void setTransientMode(bool enable, Scalar dt = 0.001, Scalar theta = 0.5);
    
    // Update old values for time stepping (called before solve)
    void updateOldValues();

private:
    // Turbulence fields
    ScalarField k;                  // Turbulent kinetic energy
    ScalarField omega;              // Specific dissipation rate
    ScalarField mu_t;               // Turbulent viscosity
    ScalarField wallDistance;       // Distance to nearest wall
    ScalarField wallShearStress;    // Wall shear stress magnitude
    
    // Gradient fields
    VectorField grad_k;             // Gradient of k
    VectorField grad_omega;         // Gradient of omega

    // Mesh references
    const std::vector<Face>& allFaces;
    const std::vector<Cell>& allCells;
    const BoundaryConditions& bcManager;
    const GradientScheme& gradientScheme;

    // Auxiliary fields for SST model
    ScalarField F1;                 // Blending function 1
    ScalarField F2;                 // Blending function 2
    ScalarField production_k;       // Production of k
    ScalarField production_omega;   // Production of omega
    ScalarField cross_diffusion;    // Cross-diffusion term

    // Model constants
    struct ModelConstants {
        // k-omega model constants (inner)
        Scalar sigma_k1 = 0.85;     // k diffusion coefficient (inner)
        Scalar sigma_omega1 = 0.5;  // omega diffusion coefficient (inner)
        Scalar beta_1 = 0.075;      // Destruction coefficient (inner)
        Scalar gamma_1 = 5.0/9.0;   // Production coefficient (inner)
        
        // k-epsilon model constants (outer) 
        Scalar sigma_k2 = 1.0;      // k diffusion coefficient (outer)
        Scalar sigma_omega2 = 0.856; // omega diffusion coefficient (outer)
        Scalar beta_2 = 0.0828;     // Destruction coefficient (outer)
        Scalar gamma_2 = 0.44;      // Production coefficient (outer)
        
        // Common constants
        Scalar C_mu = 0.09;         // Turbulent viscosity constant
        Scalar a1 = 0.31;           // Stress limiter constant
        Scalar kappa = 0.41;        // Karman constant
        Scalar E = 9.8;             // Wall function constant
        
        // Numerical constants
        Scalar beta_star = 0.09;    // Modified destruction coefficient
        Scalar sigma_d = 2.0 * sigma_omega2; // Cross-diffusion coefficient
    } constants;

    // Transient parameters
    bool enableTransient;        // Enable transient simulation
    Scalar dt;                  // Time step size
    Scalar theta;               // Crank-Nicolson parameter
    
    // Previous time step fields for transient
    ScalarField k_old;          // k from previous time step
    ScalarField omega_old;      // omega from previous time step

    // Matrix constructor for equation solving
    std::unique_ptr<Matrix> matrixConstruct;

    // Private helper methods
    void calculateBlendingFunctions(const std::vector<VectorField>& gradU,
                                   Scalar rho,
                                   Scalar mu_lam);

    /**
     * @brief Calculate production terms for k and omega
     * @param gradU Velocity gradient field
     * @param rho Fluid density
     */
    void calculateProductionTerms(const std::vector<VectorField>& gradU);

    /**
     * @brief Calculate cross-diffusion term for omega equation
     * @param rho Fluid density
     */
    void calculateCrossDiffusion(Scalar rho);

    /**
     * @brief Apply turbulence boundary conditions
     * @param fieldName Name of the turbulence field ("k" or "omega")
     * @param field Reference to the field
     * @param U_field Current velocity field
     * @param mu_lam Laminar dynamic viscosity
     * @param rho Fluid density
     */
    void applyTurbulenceBoundaryConditions(const std::string& fieldName,
                                         ScalarField& field,
                                         const VectorField& U_field,
                                         Scalar mu_lam,
                                         Scalar rho);

    /**
     * @brief Calculate y+ value for wall treatment
     * @param cellIdx Cell index
     * @param U_field Velocity field
     * @param mu_lam Laminar viscosity
     * @param rho Density
     * @return y+ value
     */
    Scalar calculateYPlus(size_t cellIdx,
                         const VectorField& U_field,
                         Scalar mu_lam,
                         Scalar rho) const;

    /**
     * @brief Limit production to prevent unrealistic values
     * @param P_k Production of k
     * @param rho Density
     * @param cellIdx Cell index
     * @return Limited production
     */
    Scalar limitProduction(Scalar P_k, Scalar rho, size_t cellIdx) const;
};

#endif