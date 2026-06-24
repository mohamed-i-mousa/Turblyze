/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TurbulenceModel.h
 * @brief Abstract base class defining the turbulence model interface
 *
 * @details TurbulenceModel is the common interface every turbulence model
 * implements. It owns no field storage; concrete models decide what state they
 * need and expose only the fields and diagnostics required by SIMPLE and
 * post-processing through this interface.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <memory>
#include <utility>
#include <vector>

// Project headers
#include "Scalar.h"
#include "StringTypes.h"
#include "Vector.h"
#include "Face.h"
#include "CellData.h"
#include "FaceData.h"

// *************************** Forward Declarations ***************************

class BoundaryConditions;
class ConvectionSchemes;
class GradientScheme;
class LinearSolver;
class Mesh;

// *************************** class TurbulenceModel **************************

class TurbulenceModel
{
public:

    using CellDataPair =
        std::vector<std::pair<NameRef, const ScalarField*>>;
    using BoundaryDataPair =
        std::vector<std::pair<NameRef, const FaceData<Scalar>*>>;
    using ResidualPair =
        std::vector<std::pair<NameRef, Scalar>>;

// ************************* Special Member Functions *************************

    /// Constructor
    TurbulenceModel() noexcept = default;

    /// Copy constructor and assignment - Not copyable
    TurbulenceModel(const TurbulenceModel&) = delete;
    TurbulenceModel& operator=(const TurbulenceModel&) = delete;

    /// Move constructor and assignment - Not movable
    TurbulenceModel(TurbulenceModel&&) = delete;
    TurbulenceModel& operator=(TurbulenceModel&&) = delete;

    /// Destructor - virtual for polymorphic deletion through the base
    virtual ~TurbulenceModel() noexcept = default;

// **************************** Runtime Selection ****************************

    /// Whether modelName selects the Laminar null-object (no turbulence)
    [[nodiscard]] static bool isLaminar(NameRef modelName) noexcept
    {
        return modelName == "Laminar";
    }

    /// Construct the turbulence model selected by name
    [[nodiscard]] static std::unique_ptr<TurbulenceModel> create
    (
        NameRef modelName,
        const Mesh& mesh,
        const BoundaryConditions& bc,
        const GradientScheme& gradScheme,
        const ConvectionSchemes& kScheme,
        LinearSolver& kSolver,
        const ConvectionSchemes& omegaScheme,
        LinearSolver& omegaSolver,
        Scalar nu,
        Scalar initialK,
        Scalar initialOmega,
        Scalar alphaK,
        Scalar alphaOmega,
        bool debug
    );

    /// Names of every selectable turbulence model
    [[nodiscard]] static NameList availableModels();

// ***************************** Turbulence Solve *****************************

    /// Solve turbulence equations for current iteration
    virtual void solve
    (
        const ScalarField&,
        const ScalarField&,
        const ScalarField&,
        const FaceFluxField&,
        const TensorField&
    ) = 0;

// ************************** Shared Accessor Methods *************************

    /// Whether the model carries turbulence (false for laminar)
    [[nodiscard]] virtual bool isTurbulent() const noexcept = 0;

    /// Get turbulent kinematic viscosity field
    [[nodiscard]] virtual const ScalarField&
    turbulentViscosity() const noexcept = 0;

    /// Get turbulent viscosity for a boundary face
    [[nodiscard]] virtual Scalar boundaryTurbulentViscosity
    (
        const Face& face,
        const BoundaryConditions&
    ) const
    {
        return turbulentViscosity()[face.ownerCell()];
    }

    /// Get kinematic wall shear stress for current state and velocity fields
    [[nodiscard]] virtual FaceData<Scalar> wallShearStress
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz
    ) const = 0;

    /// Whether wall-distance initialization converged
    [[nodiscard]] virtual bool wallDistanceConverged() const noexcept
    {
        return true;
    }

    /// Cell-centered scalar fields exported by this model (name, field)
    [[nodiscard]] virtual CellDataPair cellDataOutputs() const
    {
        return {};
    }

    /// Boundary scalar fields exported by this model (name, field)
    [[nodiscard]] virtual BoundaryDataPair boundaryDataOutputs() const
    {
        return {};
    }

    /// Residuals contributed by this model to convergence checks (name, value)
    [[nodiscard]] virtual ResidualPair residualOutputs() const
    {
        return {};
    }
};
