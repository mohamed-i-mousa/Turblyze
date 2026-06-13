/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PvdTimeSeries.h
 * @brief PVD (ParaView Data) collection file helpers for transient runs
 *
 * @details A PVD file is a small XML collection that groups per-timestep
 * `.vtu` outputs so ParaView can animate them as a single time series.
 * These helpers create and append to that XML file — no VTK library
 * dependency is involved, since the format is plain text.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Project headers
#include "Scalar.h"
#include "StringTypes.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

/// Create PVD time series file header for transient runs
void writePVDTimeSeriesHeader(const FilePath& pvdFile);

/// Append a timestep to PVD time series file
void appendPVDTimeStep
(
    const FilePath& pvdFile,
    const FilePath& vtuFile,
    Scalar timeValue
);

} // namespace VTK
