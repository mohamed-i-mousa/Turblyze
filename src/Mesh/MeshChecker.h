/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshChecker.h
 * @brief Mesh quality assessment utilities
 *
 * @details MeshChecker::check, scans every face and cell once, collecting 
 * area, volume, non-orthogonality, skewness, and aspect-ratio statistics, 
 * then prints a summary. When configured thresholds are exceeded, it give a 
 * warning.
 * Connectivity validation runs first; an invalid mesh aborts via FatalError.
 *
 * Quality metrics checked:
 * - Minimum and maximum face areas
 * - Minimum and maximum cell volumes
 * - Non-orthogonality (face-to-face angle deviation)
 * - Skewness (face/cell-centroid offset)
 * - Aspect ratio
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Mesh.h"

// *************************** namespace MeshChecker **************************

namespace MeshChecker
{

/// Run mesh quality checks and report statistics
void check(const Mesh& mesh);

}
