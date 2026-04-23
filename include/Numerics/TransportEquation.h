/******************************************************************************
 * @file TransportEquation.h
 * @brief Data bundle describing a scalar transport equation
 *
 * @details This header defines OptionalRef and TransportEquation.
 * TransportEquation bundles all data describing a scalar transport equation
 * (field, convection, diffusion, source, gradients) needed by
 * Matrix::buildMatrix().
 *
 * @struct TransportEquation
 * 
 * - Combines the field value, convection, diffusion, source,
 *   and gradient data needed by Matrix::buildMatrix().
 * 
 * - Organised by physics terms:
 *     field, transient (placeholder), convection, diffusion,
 *     source, gradient reconstruction, boundary overrides.
 *****************************************************************************/

#pragma once

#include <string>
#include <optional>

#include "CellData.h"
#include "FaceData.h"
#include "ConvectionSchemes.h"
#include "GradientScheme.h"
#include "OptionalRef.h"


struct TransportEquation
{

// Field

    /// Name of the field ("Ux", "k", "pCorr", etc.)
    std::string fieldName;

    /// Current cell-centered field values (mutable for zero-copy solve)
    ScalarField& phi;

// Transient (placeholder for future)

// Convection: div(F * phi)

    /// Face volumetric flow rates (nullopt = no convection)
    OptionalRef<FaceFluxField> flowRate = std::nullopt;

    /// Convection discretization scheme (nullopt = no convection)
    OptionalRef<ConvectionSchemes> convScheme = std::nullopt;

// Diffusion: div(Gamma * grad(phi))

    /// Cell-centered diffusion coefficient
    OptionalRef<ScalarField> Gamma = std::nullopt;

    /// Pre-interpolated face diffusion coefficient
    OptionalRef<FaceFluxField> GammaFace = std::nullopt;

// Source

    /// Explicit source term field
    const ScalarField& source;

// Gradient reconstruction

    /// Pre-computed cell gradients of phi
    const VectorField& gradPhi;

    /// Gradient reconstruction scheme
    const GradientScheme& gradScheme;

// Vector component extraction

    /// Component index for vector field equations (0=x, 1=y, 2=z)
    std::optional<int> componentIdx = std::nullopt;
};
