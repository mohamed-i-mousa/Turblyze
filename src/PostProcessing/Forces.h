/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Forces.h
 * @brief Aerodynamic force calculation on a wall patch
 *
 * @details Forces integrates the pressure and skin-friction loads over a named
 * wall patch after the SIMPLE solve converges, decomposes them onto
 * user-supplied drag and lift directions, prints the breakdown to the console,
 * and writes it to a text file alongside the VTK output. It is a read-only
 * post-processing consumer: it adds no solver state and reuses the converged
 * pressure/velocity fields and mesh geometry.
 *****************************************************************************/

#pragma once

// *************************** Forward Declarations ***************************

class SIMPLE;
class TurbulenceModel;
class Mesh;
class BoundaryConditions;
struct CaseConfiguration;

// ***************************** namespace Forces *****************************

namespace Forces
{

/// Integrate, decompose, report, and write aerodynamic forces on a wall patch
void reportForces
(
    const SIMPLE& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
);

} // namespace Forces
