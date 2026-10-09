/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MomentumTransport.h
 * @brief Abstract base class for pressure-velocity coupling algorithms
 *
 * @details MomentumTransport owns the common to every momentum transport
 * algorithm. This includes, the time control, the turbulence solve,
 * the convergence / Courant / residual evaluation, the velocity-gradient
 * reconstruction, and the solution fields. The parts that differ between
 * algorithms are reached through a virtual method.
 *
 * @class MomentumTransport
 * - Universal time control shared by steady and transient runs
 * - Runtime-selection factory (create / availableAlgorithms)
 * - Convergence and Courant-number reporting
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <memory>

// Project headers
#include "Scalar.h"
#include "StringTypes.h"
#include "Vector.h"
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "FaceData.h"
#include "GradientScheme.h"
#include "TurbulenceModel.h"

// *************************** Forward Declarations ***************************

class TimeScheme;
class ConvectionScheme;
class LinearSolver;
class Initializer;

// *************************** struct TransientFields *************************

/// Previous time step quantities for the transient runs
struct TransientFields
{
    ScalarField UxPrevStep;
    ScalarField UyPrevStep;
    ScalarField UzPrevStep;
    ScalarField UxDdtPrevStep;
    ScalarField UyDdtPrevStep;
    ScalarField UzDdtPrevStep;
    FaceFluxField fluxPrevStep;

    explicit TransientFields(const Mesh& mesh)
    :
        UxPrevStep(mesh),
        UyPrevStep(mesh),
        UzPrevStep(mesh),
        UxDdtPrevStep(mesh),
        UyDdtPrevStep(mesh),
        UzDdtPrevStep(mesh),
        fluxPrevStep(mesh)
    {}
};

// ************************* class MomentumTransport **************************

class MomentumTransport
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    MomentumTransport
    (
        const Mesh& mesh,
        BoundaryConditions& bc,
        const TimeScheme& timeScheme,
        const GradientScheme& gradScheme,
        TurbulenceModel& turbulence,
        const Initializer& initializer,
        Scalar deltaT,
        Scalar rho,
        Scalar mu,
        Count maxIterations,
        Scalar convergenceTolerance,
        Count nOuterCorrectors,
        bool debug
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    MomentumTransport(const MomentumTransport&) = delete;
    MomentumTransport& operator=(const MomentumTransport&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    MomentumTransport(MomentumTransport&&) = delete;
    MomentumTransport& operator=(MomentumTransport&&) = delete;

    /// Destructor
    virtual ~MomentumTransport() noexcept = default;

// **************************** Runtime Selection *****************************

    /// Construct the momentum transport algorithm selected by name
    [[nodiscard]] static std::unique_ptr<MomentumTransport> create
    (
        const Name& algorithm,
        const Mesh& mesh,
        BoundaryConditions& bc,
        const TimeScheme& timeScheme,
        const GradientScheme& gradScheme,
        const ConvectionScheme& momentumConvectionScheme,
        LinearSolver& momentumSolver,
        LinearSolver& pressureSolver,
        TurbulenceModel& turbulence,
        const Initializer& initializer,
        Scalar deltaT,
        Scalar rho,
        Scalar mu,
        Scalar alphaU,
        Scalar alphaP,
        Count maxIterations,
        Scalar convergenceTolerance,
        Count nNonOrthogonalCorrectors,
        Count nOuterCorrectors,
        Count nPrimeCorrectors,
        bool debug
    );

    /// Names of every selectable momentum transport algorithm
    [[nodiscard]] static NameList availableAlgorithms();

// ****************************** Solver Driver *******************************

    /// Steady: iterate to convergence. Transient: advance one time step
    void solve
    (
        Count step = 0,
        Count totalSteps = 0,
        Scalar time = S(0.0),
        TransientFields* prevStep = nullptr
    );

    /// Whether the time scheme is transient
    [[nodiscard]] bool isTransient() const noexcept;

// ***************************** Accessor Methods *****************************

    /// Get velocity field components
    [[nodiscard]] const ScalarField& Ux() const noexcept { return Ux_; }
    [[nodiscard]] const ScalarField& Uy() const noexcept { return Uy_; }
    [[nodiscard]] const ScalarField& Uz() const noexcept { return Uz_; }

    /// Get pressure field
    [[nodiscard]] const ScalarField& pressure() const noexcept { return p_; }

    /// Mesh view (nodes, faces, cells)
    [[nodiscard]] const Mesh& mesh() const noexcept
    {
        return mesh_;
    }

    /// Boundary-condition manager view
    [[nodiscard]] const BoundaryConditions& bcManager() const noexcept
    {
        return bcManager_;
    }

// ***************************** Protected Methods ****************************

protected:

    /// Run one outer-iteration
    [[nodiscard]] virtual bool outerIteration
    (
        const TransientFields* prevStep
    ) = 0;

    /// Algorithm label for banners and convergence messages
    [[nodiscard]] virtual Name algorithmName() const noexcept = 0;

    /// The convective face mass flux used by the shared driver
    [[nodiscard]] virtual const FaceFluxField&
    faceMassFlux() const noexcept
    {
        return RhieChowFlowRate_;
    }

    /// Family pressure residual for the convergence check
    [[nodiscard]] virtual Scalar pressureResidual() const noexcept = 0;

// Shared finite-volume methods

    /// Update velocity-coupled symmetry boundary coefficients
    void updateSymmetryBoundaries();

    /// Build cell and face effective viscosity
    void updateEffectiveViscosity();

    /// Interpolate the face momentum diagonal DUf_ from DU_
    void buildFaceDiagonal();

    /// Update face mass fluxes using Rhie-Chow interpolation
    void updateRhieChowFlowRate
    (
        Scalar alphaU,
        const TransientFields* prevStep = nullptr
    );

    /// Add the turbulent transpose-gradient source term: div(nuEff * (gradU)^T)
    void addTransposeGradientSource
    (
        ScalarField& UxSource,
        ScalarField& UySource,
        ScalarField& UzSource
    ) const;


// Shared driver helpers

    /// Roll the Crank-Nicolson stored velocity time derivatives forward
    void updatePrevStepDerivatives(TransientFields& prevStep);

    /// Compute limited velocity gradients into gradU_ / gradU{x,y,z}_
    void updateVelocityGradients();

    /// Solve the turbulence transport equations
    void solveTurbulence();

    /// Check convergence against the scaled-residual tolerance
    [[nodiscard]] bool checkConvergence();

    /// Compute the normalized mass-imbalance residual
    [[nodiscard]] Scalar massImbalance() const noexcept;

    /// Compute the normalized velocity residual
    [[nodiscard]] Scalar velocityResidual() const noexcept;

    /// Maximum and mean Courant number over the mesh
    struct CourantNumber
    {
        Scalar max;
        Scalar mean;
    };

    /// Compute the maximum and mean cell Courant number
    [[nodiscard]] CourantNumber computeCourant() const noexcept;

// **************************** Protected Accessors ***************************

    /// Total domain cells across every rank (cached: run-invariant)
    [[nodiscard]] Count totalDomainCells() const noexcept
    {
        return totalDomainCells_;
    }

    [[nodiscard]] BoundaryConditions& bcManager() noexcept
    {
        return bcManager_;
    }

    /// Time-derivative discretization scheme
    [[nodiscard]] const TimeScheme& timeScheme() const noexcept
    {
        return timeScheme_;
    }

    /// Gradient reconstruction scheme
    [[nodiscard]] const GradientScheme& gradientScheme() const noexcept
    {
        return gradientScheme_;
    }

    /// Turbulence model
    [[nodiscard]] const TurbulenceModel& turbulence() const noexcept
    {
        return turbulence_;
    }

    /// Kinematic viscosity
    [[nodiscard]] Scalar nu() const noexcept
    {
        return nu_;
    }

    /// Time step size [s]
    [[nodiscard]] Scalar deltaT() const noexcept
    {
        return deltaT_;
    }

    /// Whether verbose console output is enabled
    [[nodiscard]] bool debug() const noexcept
    {
        return debug_;
    }

    /// Mutable velocity field components (derived correctors write in place)
    [[nodiscard]] ScalarField& Ux() noexcept
    {
        return Ux_;
    }
    [[nodiscard]] ScalarField& Uy() noexcept
    {
        return Uy_;
    }
    [[nodiscard]] ScalarField& Uz() noexcept
    {
        return Uz_;
    }

    /// Mutable pressure field
    [[nodiscard]] ScalarField& pressure() noexcept
    {
        return p_;
    }

    /// Velocity from the previous iteration (for the velocity residual)
    [[nodiscard]] const ScalarField& UxPrevIter() const noexcept
    {
        return UxPrevIter_;
    }
    [[nodiscard]] const ScalarField& UyPrevIter() const noexcept
    {
        return UyPrevIter_;
    }
    [[nodiscard]] const ScalarField& UzPrevIter() const noexcept
    {
        return UzPrevIter_;
    }

    /// Per-component velocity gradients
    [[nodiscard]] const VectorField& gradUx() const noexcept
    {
        return gradUx_;
    }
    [[nodiscard]] const VectorField& gradUy() const noexcept
    {
        return gradUy_;
    }
    [[nodiscard]] const VectorField& gradUz() const noexcept
    {
        return gradUz_;
    }

    /// Velocity gradient tensor field
    [[nodiscard]] const TensorField& gradU() const noexcept
    {
        return gradU_;
    }

    /// Pressure gradient field
    [[nodiscard]] const VectorField& gradP() const noexcept
    {
        return gradP_;
    }
    [[nodiscard]] VectorField& gradP() noexcept
    {
        return gradP_;
    }

    /// Effective viscosity (laminar + turbulent)
    [[nodiscard]] const ScalarField& nuEff() const noexcept
    {
        return nuEff_;
    }
    [[nodiscard]] ScalarField& nuEff() noexcept
    {
        return nuEff_;
    }

    /// Effective viscosity at face centres
    [[nodiscard]] const FaceData<Scalar>& nuEffFace() const noexcept
    {
        return nuEffFace_;
    }
    [[nodiscard]] FaceData<Scalar>& nuEffFace() noexcept
    {
        return nuEffFace_;
    }

    /// Face velocity (current iteration)
    [[nodiscard]] const FaceData<Scalar>& UxAvgf() const noexcept
    {
        return UxAvgf_;
    }
    [[nodiscard]] FaceData<Scalar>& UxAvgf() noexcept
    {
        return UxAvgf_;
    }
    [[nodiscard]] const FaceData<Scalar>& UyAvgf() const noexcept
    {
        return UyAvgf_;
    }
    [[nodiscard]] FaceData<Scalar>& UyAvgf() noexcept
    {
        return UyAvgf_;
    }
    [[nodiscard]] const FaceData<Scalar>& UzAvgf() const noexcept
    {
        return UzAvgf_;
    }
    [[nodiscard]] FaceData<Scalar>& UzAvgf() noexcept
    {
        return UzAvgf_;
    }

    /// Face velocity (previous iteration)
    [[nodiscard]] const FaceData<Scalar>& UxAvgPrevIterf() const noexcept
    {
        return UxAvgPrevIterf_;
    }
    [[nodiscard]] FaceData<Scalar>& UxAvgPrevIterf() noexcept
    {
        return UxAvgPrevIterf_;
    }
    [[nodiscard]] const FaceData<Scalar>& UyAvgPrevIterf() const noexcept
    {
        return UyAvgPrevIterf_;
    }
    [[nodiscard]] FaceData<Scalar>& UyAvgPrevIterf() noexcept
    {
        return UyAvgPrevIterf_;
    }
    [[nodiscard]] const FaceData<Scalar>& UzAvgPrevIterf() const noexcept
    {
        return UzAvgPrevIterf_;
    }
    [[nodiscard]] FaceData<Scalar>& UzAvgPrevIterf() noexcept
    {
        return UzAvgPrevIterf_;
    }

    /// Mass flux through faces (Rhie-Chow)
    [[nodiscard]] const FaceFluxField& RhieChowFlowRate() const noexcept
    {
        return RhieChowFlowRate_;
    }
    [[nodiscard]] FaceFluxField& RhieChowFlowRate() noexcept
    {
        return RhieChowFlowRate_;
    }

    /// Mass flux from the previous iteration
    [[nodiscard]] const FaceFluxField& RhieChowFlowRatePrevIter() const noexcept
    {
        return RhieChowFlowRatePrevIter_;
    }
    [[nodiscard]] FaceFluxField& RhieChowFlowRatePrevIter() noexcept
    {
        return RhieChowFlowRatePrevIter_;
    }

    /// Momentum diagonal coefficients
    [[nodiscard]] const ScalarField& DU() const noexcept
    {
        return DU_;
    }
    [[nodiscard]] ScalarField& DU() noexcept
    {
        return DU_;
    }

    /// Face momentum diagonal coefficients
    [[nodiscard]] const FaceFluxField& DUf() const noexcept
    {
        return DUf_;
    }
    [[nodiscard]] FaceFluxField& DUf() noexcept
    {
        return DUf_;
    }

// ***************************** Protected Members ****************************

// Dependencies

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

    /// Reference to BCs (mutable: coefficient re-evaluation)
    BoundaryConditions& bcManager_;

    /// Time-derivative discretization scheme
    const TimeScheme& timeScheme_;

    /// Reference to gradient scheme
    const GradientScheme& gradientScheme_;

    /// Turbulence model
    TurbulenceModel& turbulence_;

// Collocated finite-volume fields

    /// Pressure gradient field
    VectorField gradP_{mesh_};

    /// Effective viscosity (laminar + turbulent)
    ScalarField nuEff_{mesh_};

    /// Effective viscosity at face centres
    FaceData<Scalar> nuEffFace_{mesh_};

    /// Face velocity (current iteration)
    FaceData<Scalar> UxAvgf_{mesh_};
    FaceData<Scalar> UyAvgf_{mesh_};
    FaceData<Scalar> UzAvgf_{mesh_};

    /// Face velocity (previous iteration)
    FaceData<Scalar> UxAvgPrevIterf_{mesh_};
    FaceData<Scalar> UyAvgPrevIterf_{mesh_};
    FaceData<Scalar> UzAvgPrevIterf_{mesh_};

    /// Mass flux through faces (Rhie-Chow)
    FaceFluxField RhieChowFlowRate_{mesh_};

    /// Mass flux from the previous iteration
    FaceFluxField RhieChowFlowRatePrevIter_{mesh_};

    /// Momentum diagonal coefficients
    ScalarField DU_{mesh_};

    /// Face momentum diagonal coefficients
    FaceFluxField DUf_{mesh_};

// ****************************** Private Members *****************************

private:

// Physical properties

    /// Kinematic viscosity
    Scalar nu_;

// Algorithm parameters

    /// Time step size [s]
    Scalar deltaT_;

    /// Maximum number of steady outer iterations
    Count maxIterations_;

    /// Fixed number of outer correctors per transient time step
    Count nOuterCorrectors_;

    /// Convergence tolerance
    Scalar tolerance_;

    /// Enable verbose console output
    bool debug_;

    /// Whether the inner loop should print per-iteration residual lines
    bool reportPerIteration_ = true;

// Solution fields

    /// Velocity fields
    ScalarField Ux_{mesh_};
    ScalarField Uy_{mesh_};
    ScalarField Uz_{mesh_};

    /// Pressure field
    ScalarField p_{mesh_};

    /// Velocity from previous iteration (for the velocity residual)
    ScalarField UxPrevIter_{mesh_};
    ScalarField UyPrevIter_{mesh_};
    ScalarField UzPrevIter_{mesh_};

// Gradient fields

    /// Per-component velocity gradients
    VectorField gradUx_{mesh_};
    VectorField gradUy_{mesh_};
    VectorField gradUz_{mesh_};

    /// Velocity gradient tensor field
    TensorField gradU_{mesh_};

    /// Total domain cells across every rank (reduced once at construction)
    Count totalDomainCells_ = 0;

// Residual tracking for convergence

    /// First-iteration reference values for scaled residuals
    Scalar massImbalance0_ = S(0.0);
    Scalar velocityResidual0_ = S(0.0);
    Scalar pressureResidual0_ = S(0.0);
    ScalarList turbulenceResidual0_;

    /// Most recent scaled residuals
    Scalar lastScaledMass_ = S(0.0);
    Scalar lastScaledVelocity_ = S(0.0);
    Scalar lastScaledPressure_ = S(0.0);
    TurbulenceModel::ResidualPair lastScaledTurbulence_;
};