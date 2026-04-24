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

#include <string>

#include "Scalar.h"


namespace VTK
{

/**
 * @brief Create PVD time series file header for transient runs
 * @param pvdFilename Path to .pvd file to create
 */
void writePVDTimeSeriesHeader(const std::string& pvdFilename);

/**
 * @brief Append a timestep to PVD time series file
 * @param pvdFilename Path to existing .pvd file
 * @param vtuFilename Relative path to .vtu file for this timestep
 * @param timeValue Physical time value for this timestep
 */
void appendPVDTimeStep
(
    const std::string& pvdFilename,
    const std::string& vtuFilename,
    Scalar timeValue
);

} // namespace VTK
