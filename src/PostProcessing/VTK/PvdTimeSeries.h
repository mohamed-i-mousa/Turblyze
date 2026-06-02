/******************************************************************************
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

// Standard library headers
#include <string>

// Project headers
#include "Scalar.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

/// Create PVD time series file header for transient runs
void writePVDTimeSeriesHeader(const std::string& pvdFilename);

/// Append a timestep to PVD time series file
void appendPVDTimeStep
(
    const std::string& pvdFilename,
    const std::string& vtuFilename,
    Scalar timeValue
);

} // namespace VTK
