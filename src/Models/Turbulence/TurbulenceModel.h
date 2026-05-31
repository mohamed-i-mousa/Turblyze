/******************************************************************************
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
#include <string_view>
#include <utility>
#include <vector>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "CellData.h"
#include "FaceData.h"

// *************************** Forward Declarations ***************************

class BoundaryConditions;

// *************************** class TurbulenceModel **************************

class TurbulenceModel
{
public:

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

    /// Whether wall-distance initialization converged
    [[nodiscard]] virtual bool wallDistanceConverged() const noexcept
    {
        return true;
    }

    /// Cell-centered scalar fields exported by this model (name, field)
    [[nodiscard]] virtual
    std::vector<std::pair<std::string_view, const ScalarField*>>
    cellDataOutputs() const
    {
        return {};
    }

    /// Boundary scalar fields exported by this model (name, field)
    [[nodiscard]] virtual
    std::vector<std::pair<std::string_view, const FaceData<Scalar>*>>
    boundaryDataOutputs() const
    {
        return {};
    }

    /// Residuals contributed by this model to convergence checks (name, value)
    [[nodiscard]] virtual std::vector<std::pair<std::string_view, Scalar>>
    residualOutputs() const
    {
        return {};
    }
};
