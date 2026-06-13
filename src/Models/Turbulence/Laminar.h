/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Laminar.h
 * @brief Laminar (no-turbulence) null-object turbulence model
 *
 * @details Laminar satisfies the TurbulenceModel interface for runs without a
 * turbulence model. Its turbulent viscosity is zero, solve() is a no-op, and
 * isTurbulent() reports false so the solver skips turbulence-specific source
 * terms and reporting. It computes laminar wall shear stress on demand and
 * never computes wall distance.
 *
 * @class Laminar
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <algorithm>
#include <cmath>

// Project headers
#include "Mesh.h"
#include "TurbulenceModel.h"

// ******************************* class Laminar ******************************

class Laminar final : public TurbulenceModel
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    Laminar(const Mesh& mesh, Scalar nu) noexcept
    :
        mesh_{mesh},
        nu_{nu}
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

    /// No turbulence equations are solved for laminar runs
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

    /// Compute kinematic wall shear stress magnitude (tau/rho) on wall faces
    [[nodiscard]] FaceData<Scalar> wallShearStress
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz
    ) const override
    {
        FaceData<Scalar> shearStress(S(0.0));

        for (const Face& face : mesh_.faces())
        {
            if (!face.isBoundary())
            {
                continue;
            }

            const auto& patch = face.patch();
            if
            (
                !patch.has_value()
             || patch->get().type() != PatchType::wall
            )
            {
                continue;
            }

            const Index cellIdx = face.ownerCell();
            const Vector& normal = face.normal();
            const Vector cellVelocity(Ux[cellIdx], Uy[cellIdx], Uz[cellIdx]);
            const Scalar normalVelocity = dot(cellVelocity, normal);
            const Vector tangentVelocity =
                cellVelocity - normalVelocity * normal;
            const Scalar wallDistance =
                std::max
                (
                    std::abs(dot(face.dPf(), normal)),
                    vSmallValue
                );

            shearStress[face.idx()] =
                nu_ * magnitude(tangentVelocity) / wallDistance;
        }

        return shearStress;
    }

// ****************************** Private Members *****************************

private:

    /// Mesh view for laminar wall-shear evaluation
    const Mesh& mesh_;

    /// Laminar kinematic viscosity
    Scalar nu_;

    /// Zero turbulent viscosity field
    ScalarField nut_{S(0.0)};

};
