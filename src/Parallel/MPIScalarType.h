/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MPIScalarType.h
 * @brief The MPI datatypes matching project primitive types
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// External library headers
#include <mpi.h>

// Project headers
#include "Scalar.h"

// ****************************** MPI Data Types ******************************

/// MPI datatype of one Scalar
[[nodiscard]] inline MPI_Datatype MPIScalarType()
{
    return sizeof(Scalar) == 8 ? MPI_DOUBLE : MPI_FLOAT;
}
