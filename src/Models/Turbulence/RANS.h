/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file RANS.h
 * @brief Abstract base for two-equation Reynolds-Averaged Navier-Stokes models
 *
 * @details RANS sits between the root TurbulenceModel interface and concrete
 * two-equation models. It owns the state and services common to such models:
 * mesh/BC views, equation assembly services, turbulent viscosity, k transport
 * state, wall-function geometry/diagnostics, and residual bookkeeping.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <memory>
#include <utility>
#include <vector>

// Project headers
#include "Scalar.h"
#include "Integer.h"
#include "StringTypes.h"
#include "Vector.h"
#include "CellData.h"
#include "FaceData.h"
#include "TurbulenceModel.h"

// *************************** Forward Declarations ***************************

class Mesh;
class BoundaryConditions;
class GradientScheme;
class ConvectionSchemes;
class LinearSolver;
class Matrix;
class Face;
enum class BCType;
enum class Field;

// ******************************** class RANS ********************************

class RANS : public TurbulenceModel
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    RANS
    (
        const Mesh& mesh,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionSchemes& kScheme,
        LinearSolver& kSolver,
        const ConvectionSchemes& dissipationScheme,
        LinearSolver& dissipationSolver,
        Scalar nu,
        Scalar alphaK,
        Scalar alphaDissipation,
        bool debug
    );

    /// Copy constructor and assignment - Not copyable (reference members)
    RANS(const RANS&) = delete;
    RANS& operator=(const RANS&) = delete;

    /// Move constructor and assignment - Not movable (reference members)
    RANS(RANS&&) = delete;
    RANS& operator=(RANS&&) = delete;

    /// Destructor - out-of-line to anchor the unique_ptr<Matrix> member
    ~RANS() noexcept override;

// ***************************** Turbulence Solve *****************************

    /// Solve turbulence equations for current iteration
    void solve
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz,
        const FaceFluxField& flowRateFace,
        const TensorField& gradU
    ) override = 0;

// ************************ Inlet Condition Calculators ***********************

    /// Calculate inlet/initial turbulent kinetic energy
    [[nodiscard]] static Scalar inletK
    (
        const Vector& velocity,
        Scalar turbulenceIntensity
    ) noexcept;

// ***************************** Accessor Methods *****************************

    /// Get turbulent kinematic viscosity field
    [[nodiscard]] const ScalarField&
    turbulentViscosity() const noexcept override
    {
        return nut_;
    }

    /// Get turbulent viscosity for a boundary face
    [[nodiscard]] Scalar boundaryTurbulentViscosity
    (
        const Face& face,
        const BoundaryConditions& bcManager
    ) const override;

    /// Wall shear stress (tau/rho) from the wall-function state and velocity
    [[nodiscard]] FaceData<Scalar> wallShearStress
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz
    ) const override;

    /// Get turbulent kinetic energy field
    [[nodiscard]] const ScalarField& k() const noexcept
    {
        return k_;
    }

    /// Get the model's dissipation field
    [[nodiscard]] virtual const ScalarField& dissipation() const noexcept = 0;

    /// Name of the dissipation field for output labelling ("omega"/"epsilon")
    [[nodiscard]] virtual NameRef dissipationName() const noexcept = 0;

    /// Get normalised k change from the most recent solve
    [[nodiscard]] Scalar lastKResidual() const noexcept
    {
        return lastKResidual_;
    }

    /// Get normalised dissipation change from the most recent solve
    [[nodiscard]] Scalar lastDissipationResidual() const noexcept
    {
        return lastDissipationResidual_;
    }

    /// Whether wall-distance initialization converged
    [[nodiscard]] bool wallDistanceConverged() const noexcept override
    {
        return wallDistanceConverged_;
    }

    /// Cell-centered scalar fields exported by RANS models (name, field)
    [[nodiscard]] CellDataPair cellDataOutputs() const override;

    /// Boundary scalar fields exported by RANS models (name, field)
    [[nodiscard]] BoundaryDataPair boundaryDataOutputs() const override;

    /// Residuals contributed by RANS models to convergence checks (name, value)
    [[nodiscard]] ResidualPair residualOutputs() const override;

// ****************************** Shared Members ******************************

protected:

// Dependencies and services

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

    /// Reference to BCs
    const BoundaryConditions& bcManager_;

    /// Reference to gradient scheme
    const GradientScheme& gradientScheme_;

    /// Matrix constructor
    std::unique_ptr<Matrix> matrixConstruct_;

    /// Reference to k convection scheme
    const ConvectionSchemes& kConvectionScheme_;

    /// Linear solver for k equation
    LinearSolver& kSolver_;

    /// Reference to dissipation convection scheme
    const ConvectionSchemes& dissipationConvectionScheme_;

    /// Linear solver for dissipation equation
    LinearSolver& dissipationSolver_;

// Physical and algorithm parameters

    /// Laminar kinematic viscosity
    Scalar nu_;

    /// Under-relaxation factor for k equation
    Scalar alphaK_;

    /// Under-relaxation factor for dissipation equation
    Scalar alphaDissipation_;

    /// Enable verbose console output
    bool debug_;

// Common transport state

    /// Turbulent kinematic viscosity
    ScalarField nut_{S(0.0)};

    /// Turbulent kinetic energy
    ScalarField k_{S(1e-6)};

    /// Previous-iteration k snapshot for residual computation
    ScalarField kPrev_;

    /// Cell gradient of k
    VectorField gradK_;

    /// Normalised k change from the most recent solve
    Scalar lastKResidual_ = S(1e9);

    /// Normalised dissipation change from the most recent solve
    Scalar lastDissipationResidual_ = S(1e9);

// Wall distance and wall-function state

    /// Distance to nearest wall
    ScalarField wallDistance_{S(1.0)};

    /// Coordinates of the nearest wall point per cell (for mesh-wave)
    VectorField nearestWallPoint_;

    /// meshWave wall-distance loop convergence flag
    bool wallDistanceConverged_ = false;

    /// Owner-cell to wall-face perpendicular distance
    FaceData<Scalar> y_;

    /// y+
    FaceData<Scalar> yPlus_;

    /// Wall-function nut values on wall faces
    FaceData<Scalar> nutWall_;

    /// Area-based weight per wall face (face area / total wall area of cell)
    FaceData<Scalar> wallFaceWeight_;

    /// Indices into mesh_.faces for faces with model wall-function BCs
    IndexList wallFunctionFaceIndices_;

    /// Unique cell indices adjacent to wall-function faces
    IndexList wallCellIndices_;

    /// Wall-to-total boundary area fraction per wall cell
    ScalarList wallCellFraction_;

    /// y+ crossover between viscous sublayer and log region
    Scalar yPlusLam_ = S(11.225);

// ****************************** Shared Methods ******************************

    /// Map cell-centered diffusion coefficients to faces for assembly
    void cellToFaceDiffusion
    (
        const ScalarField& cellGamma,
        FaceFluxField& faceGamma
    ) const;

    /// Update wall distance field using mesh-wave coordinate propagation
    void updateWallDistance();

    /// Compute y+ crossover via fixed-point iteration
    void updateYPlusLam(Scalar kappa, Scalar E);

    /// Build wall-function face lists and area weights
    void initializeWallFunctionGeometry
    (
        const BoundaryConditions& bcManager,
        Field wallFunctionField,
        BCType wallFunctionType
    );

    /// Model-specific Cμ^0.25 used by the wall functions
    [[nodiscard]] virtual Scalar cmu25() const noexcept = 0;

    /// Update y+ field on wall-function faces
    void updateYPlus();

    /// Compute strain-rate magnitude: ||S|| = sqrt(2 S_ij S_ij)
    [[nodiscard]] ScalarField computeStrainRateMagnitude
    (
        const TensorField& gradU
    ) const;

    /// Compute cell velocity divergence, allocated on demand
    [[nodiscard]] ScalarField velocityDivergence
    (
        const FaceFluxField& flowRateFace
    ) const;

    /// Update both residuals against the pre-solve snapshots
    void updateResiduals
    (
        const ScalarField& dissipation,
        const ScalarField& dissipationPrev
    );

    /// Compute normalised field change against a previous snapshot
    [[nodiscard]] Scalar normalisedFieldResidual
    (
        const ScalarField& field,
        const ScalarField& previousField
    ) const;
};
