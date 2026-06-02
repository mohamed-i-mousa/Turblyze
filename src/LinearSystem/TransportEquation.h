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

// ********************************** Headers *********************************

// Standard library headers
#include <optional>

// Project headers
#include "CellData.h"
#include "Field.h"
#include "FaceData.h"
#include "ConvectionSchemes.h"
#include "GradientScheme.h"
#include "OptionalRef.h"

// ************************* struct TransportEquation *************************

struct TransportEquation
{

// *********************************** Field **********************************

    /// Field identifier
    Field field;

    /// Current cell-centered field values (mutable for zero-copy solve)
    ScalarField& phi;

// ********************************* Transient ********************************

// (placeholder for future)

// ******************************** Convection ********************************

// Convection: div(F * phi)

    /// Face volumetric flow rates (nullopt = no convection)
    OptionalRef<FaceFluxField> flowRate = std::nullopt;

    /// Convection discretization scheme (nullopt = no convection)
    OptionalRef<ConvectionSchemes> convScheme = std::nullopt;

// ********************************* Diffusion ********************************

// Diffusion: div(Gamma * grad(phi))

    /// Cell-centered diffusion coefficient
    OptionalRef<ScalarField> Gamma = std::nullopt;

    /// Pre-interpolated face diffusion coefficient
    OptionalRef<FaceFluxField> GammaFace = std::nullopt;

// ********************************** Source **********************************

    /// Explicit source term field
    const ScalarField& source;

// ************************** Gradient Reconstruction *************************

    /// Pre-computed cell gradients of phi
    const VectorField& gradPhi;

    /// Gradient reconstruction scheme
    const GradientScheme& gradScheme;
};
