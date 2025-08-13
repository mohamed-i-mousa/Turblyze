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


class SIMPLE {
public:
    SIMPLE(const std::vector<Face>& faces,
           const std::vector<Cell>& cells, 
           const BoundaryConditions& bc,
           const GradientScheme& gradScheme,
           const ConvectionScheme& convScheme);

    // Main SIMPLE algorithm steps
    void solve();
    
    // Individual SIMPLE steps
    void solveMomentumEquations();
    void calculateFaceFluxes();
    void solvePressureCorrection();
    void correctVelocity();
    void correctPressure();
    void correctMassFluxes();
    bool checkConvergence();
    bool checkDivergence();

    // Rhie-Chow interpolation
    void calculateRhieChowFaceVelocities();
    void calculateMassFluxes();

    // Getters for solution fields
    const VectorField& getVelocity() const { return U; }
    const ScalarField& getPressure() const { return p; }
    const FaceFluxField& getMassFlux() const { return massFlux; }
    
    // 3D-specific getters
    ScalarField getVelocityX() const;
    ScalarField getVelocityY() const;
    ScalarField getVelocityZ() const;

    // Setters for algorithm parameters
    void setRelaxationFactors(Scalar alpha_U, Scalar alpha_p);
    void setConvergenceTolerance(Scalar tol);
    void setMaxIterations(int maxIter);
    void enableTurbulenceModeling(bool enable = true);
    void setPhysicalProperties(Scalar rho_in, Scalar mu_in) { rho = rho_in; mu = mu_in; }
    
    // Turbulence getters
    const ScalarField* getTurbulentKineticEnergy() const;
    const ScalarField* getSpecificDissipationRate() const;
    const ScalarField* getTurbulentViscosity() const;
    const ScalarField* getWallDistance() const;

private:
    // Mesh references
    const std::vector<Face>& allFaces;
    const std::vector<Cell>& allCells;
    const BoundaryConditions& bcManager;
    const GradientScheme& gradientScheme;
    const ConvectionScheme& convectionScheme;

    // Physical properties
    Scalar rho;           // Density
    Scalar mu;            // Dynamic viscosity

    // Algorithm parameters
    Scalar alpha_U;       // Under-relaxation factor for velocity
    Scalar alpha_p;       // Under-relaxation factor for pressure
    int maxIterations;    // Maximum number of iterations
    Scalar tolerance;     // Convergence tolerance
    bool enableTurbulence;        // Enable turbulence modeling
    
    // Turbulence model
    std::unique_ptr<KOmegaSST> turbulenceModel;

    // Solution fields
    VectorField U;        // Velocity field
    ScalarField p;        // Pressure field
    ScalarField p_prime;  // Pressure correction field
    Scalar lastPressureCorrectionRMS = S(1e9); // track p' RMS before reset
    
    // Face-based fields for Rhie-Chow interpolation
    FaceVectorField U_face;      // Face velocities
    FaceFluxField massFlux;      // Mass flux through faces
    FaceFluxField volumeFlux;    // Volume flux through faces

     // Previous-iteration fields (for under-relaxation effects at faces)
     VectorField U_prev;          // Cell-centered velocity from previous iteration
     FaceVectorField U_face_prev; // Face velocity from previous iteration

    // Momentum equation coefficients per component (needed for Rhie-Chow)
    ScalarField a_Ux;     // Diagonal coefficients of U_x momentum equation
    ScalarField a_Uy;     // Diagonal coefficients of U_y momentum equation
    ScalarField a_Uz;     // Diagonal coefficients of U_z momentum equation
    VectorField H_U;      // H/A terms for momentum equations
    
    // Gradient fields
    VectorField gradP;    // Pressure gradient
    
    // Matrix constructor and solver objects
    std::unique_ptr<Matrix> matrixConstruct;
    
    // Helper functions
    void initialize(const Vector& initialVelocity, Scalar initialPressure);
    
    // Utilities
    Scalar calculateMassImbalance() const;
    Scalar calculateVelocityResidual() const;
    Scalar calculatePressureResidual() const;
     void printSolutionStatistics();

    // Linear interpolation utilities used in Rhie-Chow and elsewhere
    void computeLinearWeights(const Face& face, Scalar& w_P, Scalar& w_N) const;
    Vector linearInterpolation(const Face& face, const VectorField& cellField) const;
    Scalar linearInterpolation(const Face& face, const ScalarField& cellField) const;
};

#endif