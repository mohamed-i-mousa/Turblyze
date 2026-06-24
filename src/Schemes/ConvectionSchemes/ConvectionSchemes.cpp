/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ConvectionSchemes.cpp
 * @brief Runtime selection of convection schemes
 *****************************************************************************/

// ********************************** Headers *********************************

#include "ConvectionSchemes.h"

// Project headers
#include "CentralDifferenceScheme.h"
#include "RuntimeSelection.h"
#include "SecondOrderUpwindScheme.h"
#include "UpwindScheme.h"

// **************************** Runtime Selection ****************************

std::unique_ptr<ConvectionSchemes> ConvectionSchemes::create(Name schemeName)
{
    if (schemeName == "Upwind")
    {
        return std::make_unique<UpwindScheme>();
    }

    if (schemeName == "CentralDifference")
    {
        return std::make_unique<CentralDifferenceScheme>();
    }

    if (schemeName == "SecondOrderUpwind")
    {
        return std::make_unique<SecondOrderUpwindScheme>();
    }

    RuntimeSelection::unknownSelection
    (
        "convection scheme",
        schemeName,
        availableSchemes()
    );
}


NameList ConvectionSchemes::availableSchemes()
{
    return {"Upwind", "CentralDifference", "SecondOrderUpwind"};
}
