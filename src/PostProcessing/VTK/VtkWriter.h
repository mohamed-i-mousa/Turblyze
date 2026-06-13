/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file VtkWriter.h
 * @brief VTK UnstructuredGrid (.vtu) writer for 3D volume cell data
 *
 * @details Exports the 3D unstructured mesh as VTK polyhedron cells together
 * with cell-centered scalar and vector fields in the VTK UnstructuredGrid
 * format, suitable for ParaView volume rendering, slicing, clipping, and
 * isosurface extraction.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <array>
#include <map>

// Project headers
#include "Scalar.h"
#include "Mesh.h"
#include "CellData.h"
#include "StringTypes.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

/// Alias for cell data maps
using ScalarFieldMap = std::map<Name, const ScalarField*>;
using VectorFieldMap = std::map<Name, std::array<const ScalarField*, 3>>;

/// Write simulation results to VTK UnstructuredGrid (.vtu) file
void writeVtkUnstructuredGrid
(
    const FilePath& filename,
    const Mesh& mesh,
    const ScalarFieldMap& scalarFields = {},
    const VectorFieldMap& vectorFields = {},
    bool debug = false
);

} // namespace VTK
