/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Segregated.h
 * @brief Abstract base class for segregated pressure-correction algorithms
 *
 * @details Segregated owns the pressure-correction and Rhie-Chow
 * shared by every segregated coupling algorithm (SIMPLE, PISO):
 * implicit momentum assembly, the pressure-correction Poisson solve, the 
 * velocity and face-flux corrections, and the momentum diagonal coefficients.
 * Concrete algorithms derive from it and implement only
 * outerIteration() and algorithmName()
 *
 * @class Segregated
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <optional>

// Project headers
#include "MomentumTransport.h"
#include "ConvectionScheme.h"
#include "Matrix.h"
#include "LinearSolvers.h"
#include "TransportEquation.h"

// *************************** Forward Declarations ***************************

class Initializer;

// ***************************** class Segregated *****************************

class Segregated : public MomentumTransport
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    Segregated
    (
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
        bool debug
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    Segregated(const Segregated&) = delete;
    Segregated& operator=(const Segregated&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    Segregated(Segregated&&) = delete;
    Segregated& operator=(Segregated&&) = delete;

    /// Destructor
    ~Segregated() noexcept override = default;

// ***************************** Protected Members ****************************

protected:

    /// Pressure correction RMS normalized by the pressure RMS
    [[nodiscard]] Scalar pressureResidual() const noexcept override;

    /// Build the gradient sources and reset the DU_ accumulators
    void assembleMomentum();

    /// Solve the three momentum components
    void solveMomentum(const TransientFields* prevStep);

    /// Compute the momentum diagonal coefficient DU_ = 3 V / sum(a_P)
    void diagonalDU(Index momentumComponent);

    /// Update face mass fluxes using Rhie-Chow interpolation
    void updateRhieChowFlowRate(const TransientFields* prevStep)
    {
        MomentumTransport::updateRhieChowFlowRate(alphaU_, prevStep);
    }

    /// Assemble and solve the pressure correction equation
    void solvePressureCorrection();

    /// Apply velocity correction: U = U* - D*gradPCorr
    void correctVelocity();

    /// Update face mass fluxes from the pressure correction gradient
    void correctFlowRate();

    /// Update pressure with under-relaxation and reset pCorr
    void correctPressure();

    /// Add the transpose-gradient source to the momentum source terms
    void addTransposeGradientSource()
    {
        MomentumTransport::addTransposeGradientSource
        (
            UxSource_,
            UySource_,
            UzSource_
        );
    }


    /// Build the transient term for one velocity component
    [[nodiscard]] std::optional<TransientTerm> ddtTerm
    (
        const TransientFields* prevStep,
        ScalarField TransientFields::* phiPrevStep,
        ScalarField TransientFields::* ddtPrevStep
    ) const;

// **************************** Protected Accessors ***************************

    /// Momentum convection scheme
    [[nodiscard]] const ConvectionScheme&
    momentumConvectionScheme() const noexcept
    {
        return momentumConvectionScheme_;
    }

    /// Matrix constructor and solver object
    [[nodiscard]] const Matrix& matrixConstruct() const noexcept
    {
        return matrixConstruct_;
    }
    [[nodiscard]] Matrix& matrixConstruct() noexcept
    {
        return matrixConstruct_;
    }


    /// Momentum source terms
    [[nodiscard]] const ScalarField& UxSource() const noexcept
    {
        return UxSource_;
    }
    [[nodiscard]] const ScalarField& UySource() const noexcept
    {
        return UySource_;
    }
    [[nodiscard]] const ScalarField& UzSource() const noexcept
    {
        return UzSource_;
    }

// ****************************** Private Members *****************************

private:

// Segregated dependencies

    /// Reference to the momentum convection scheme
    const ConvectionScheme& momentumConvectionScheme_;

    /// Linear solver for momentum equations
    LinearSolver& momentumSolver_;

    /// Linear solver for the pressure correction equation
    LinearSolver& pressureSolver_;

    /// Matrix constructor and solver object
    Matrix matrixConstruct_;

// Algorithm parameters

    /// Under-relaxation factor for velocity
    Scalar alphaU_;

    /// Under-relaxation factor for pressure
    Scalar alphaP_;

    /// Non-orthogonal corrector sub-iterations for the p' equation
    Count nNonOrthogonalCorrectors_;

// Pressure-correction fields

    /// Pressure correction field
    ScalarField pCorr_{mesh_};

    /// Pressure correction gradient field
    VectorField gradPCorr_{mesh_};

    /// True when p' has no Dirichlet anchor on any rank (pure Neumann)
    bool pCorrNeedsNullSpace_ = false;

// Momentum assembly fields


    /// Momentum source terms
    ScalarField UxSource_{mesh_};
    ScalarField UySource_{mesh_};
    ScalarField UzSource_{mesh_};

// Pressure-correction assembly fields

    /// Mass imbalance source for the pressure correction equation
    ScalarField massImbalanceSrc_{mesh_};

    /// Track pressure correction RMS before reset
    Scalar lastPressureCorrectionRMS_ = S(1e9);
};