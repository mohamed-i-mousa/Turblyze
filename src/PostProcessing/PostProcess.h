/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PostProcess.h
 * @brief After-solve reporting and result export
 *****************************************************************************/

#pragma once

// *************************** Forward Declarations ***************************

class TurbulenceModel;
class Mesh;
class SIMPLE;
struct CaseConfiguration;

// *************************** namespace PostProcess **************************

namespace PostProcess
{

/// Extract solution fields and print flow statistics
void reportStatistics(const SIMPLE& solver);

/// Write VTK volume and wall-boundary results
void exportResults
(
    const SIMPLE& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const CaseConfiguration& config
);

} // namespace PostProcess
