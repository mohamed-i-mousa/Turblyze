/******************************************************************************
 * @file Laminar.h
 * @brief Laminar (no-turbulence) null-object turbulence model
 *
 * @details Laminar satisfies the TurbulenceModel interface for runs without a
 * turbulence model. Its turbulent viscosity is zero, solve() is a no-op, and
 * isTurbulent() reports false so the solver skips turbulence-specific source
 * terms and reporting. It never computes wall distance.
 *
 * @class Laminar
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "TurbulenceModel.h"

// *************************** Forward Declarations ***************************

class Mesh;

// ******************************* class Laminar ******************************

class Laminar final : public TurbulenceModel
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    Laminar(const Mesh&, Scalar)
    {}

    /// Copy constructor and assignment - Not copyable (non-copyable base)
    Laminar(const Laminar&) = delete;
    Laminar& operator=(const Laminar&) = delete;

    /// Move constructor and assignment - Not movable (non-copyable base)
    Laminar(Laminar&&) = delete;
    Laminar& operator=(Laminar&&) = delete;

    /// Destructor
    ~Laminar() noexcept override = default;

// ***************************** Turbulence Solve *****************************

    /// Laminar model has no turbulence equations to solve
    void solve
    (
        const ScalarField&,
        const ScalarField&,
        const ScalarField&,
        const FaceFluxField&,
        const TensorField&
    ) override
    {}

// ***************************** Accessor Methods *****************************

    /// Turbulent kinematic viscosity is zero for laminar runs
    [[nodiscard]] const ScalarField&
    turbulentViscosity() const noexcept override
    {
        return nut_;
    }

    /// A laminar model carries no turbulence
    [[nodiscard]] bool isTurbulent() const noexcept override
    {
        return false;
    }

// ****************************** Private Members *****************************

private:

    /// Zero turbulent viscosity field
    ScalarField nut_{S(0.0)};
};
