/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Vector.cpp
 * @brief Stream output for the Vector class
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Vector.h"

// Standard library headers
#include <ostream>
#include <iomanip>

// *************************** Non-Member Functions ***************************

std::ostream& operator<<(std::ostream& os, const Vector& p)
{
    // save the current format 
    const auto flags = os.flags();
    const auto prec = os.precision();

    // change the format for vector output
    os  << std::fixed << std::setprecision(6);

    // vector output
    os  << '(' << p.x() << ", " << p.y() << ", " << p.z() << ')';

    // restore the original format
    os.flags(flags);
    os.precision(prec);

    return os;
}
