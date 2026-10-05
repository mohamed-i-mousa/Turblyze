/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BCLoader.cpp
 * @brief Case-file boundary condition registration
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "BCLoader.h"

// Standard library headers
#include <format>
#include <iostream>
#include <map>
#include <memory>

// Project headers
#include "BoundaryTypeFactory.h"
#include "CaseReader.h"
#include "ErrorHandler.h"
#include "FixedValue.h"
#include "Logger.h"
#include "Reduce.h"
#include "RuntimeSelection.h"
#include "TurbulenceModel.h"
#include "Vector.h"
#include "ZeroGradient.h"
#include "kOmegaSST.h"

// **************************** namespace BCLoader ****************************

namespace BCLoader
{

// ***************************** Internal Helpers *****************************

namespace
{

[[noreturn]] void unknownTypeToken
(
    const Name& bcType,
    const Name& fieldName,
    const Name& patchName,
    const Message& validList
)
{
    FatalError
    (
        "Unknown boundary condition type '" + bcType
      + "' for field '" + fieldName
      + "' on patch '" + patchName
      + "'. Valid types: " + validList
    );
}


/// Validate the type token against the field's selectable names, then
/// construct and register the boundary condition
void registerBC
(
    BoundaryConditions::BCs& bcs,
    const Name& patchName,
    Field field,
    const Name& bcType,
    Scalar value,
    const Name& sectionName
)
{
    if
    (
        !RuntimeSelection::isKnown
        (
            bcType,
            BoundaryTypeFactory::availableTypes(field)
        )
    )
    {
        unknownTypeToken
        (
            bcType,
            sectionName,
            patchName,
            RuntimeSelection::joinNames
            (
                BoundaryTypeFactory::availableTypes(field)
            )
        );
    }

    bcs[patchName][field] =
        BoundaryTypeFactory::create(bcType, field, value);
}


void validateWallFunctionSetup
(
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    if (!TurbulenceModel::isRANS(config.turbulenceModel))
    {
        return;
    }

    for (const auto& patch : mesh.patches())
    {
        if (patch.type() == PatchType::processor
         || !bcManager.boundaryType(patch.name(), Field::Ux).isWall())
        {
            continue;
        }

        const Name& patchName = patch.name();

        const bool kIsWF =
            bcManager.boundaryType(patchName, Field::k).isWallModelled();
        const bool omegaIsWF =
            bcManager.boundaryType(patchName, Field::omega).isWallModelled();
        const bool nutIsWF =
            bcManager.boundaryType(patchName, Field::nut).isWallModelled();

        const int wfCount = int(kIsWF) + int(omegaIsWF) + int(nutIsWF);

        if (wfCount == 0 || wfCount == 3)
        {
            continue;
        }

        FatalError
        (
            "Wall patch '" + patchName
          + "': wall functions must be configured as a complete triplet "
            "(k + omega + nut) or omitted entirely. Found: k="
          + (kIsWF     ? "WF" : "non-WF")
          + ", omega=" + (omegaIsWF ? "WF" : "non-WF")
          + ", nut="   + (nutIsWF   ? "WF" : "non-WF") + "."
        );
    }
}

} // namespace


// *********************************** Load ***********************************

BoundaryConditions load
(
    const CaseReader& reader,
    const CaseConfiguration& config,
    const Mesh& mesh
)
{
    std::cout << '\n';
    Logger::sectionHeader("Setting Boundary Conditions");

    BoundaryConditions::BCs bcs;

    for (const auto& face : mesh.faces())
    {
        if (face.isBoundary() && face.patch() == nullptr)
        {
            FatalError
            (
                std::format
                (
                    "Boundary face {} has no patch after linking.",
                    face.idx()
                )
            );
        }
    }

    const auto& BCs = reader.section("boundaryConditions");

    if (BCs.hasSection("U"))
    {
        const auto& velocityBCs = BCs.section("U");

        for (const auto& patchName : velocityBCs.sectionNames())
        {
            const auto& patchBC = velocityBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            Vector value{};
            if (bcType == "fixedValue")
            {
                value = patchBC.lookup<Vector>("value");
            }
            else if (bcType == "fixedGradient")
            {
                value = patchBC.lookup<Vector>("gradient");
            }

            // The case-file vector entry fans out into the scalar components
            registerBC(bcs, patchName, Field::Ux, bcType, value.x(), "U");
            registerBC(bcs, patchName, Field::Uy, bcType, value.y(), "U");
            registerBC(bcs, patchName, Field::Uz, bcType, value.z(), "U");
        }
    }

    bool hasFixedPressure = false;

    if (BCs.hasSection("p"))
    {
        const auto& pressureBCs = BCs.section("p");

        for (const auto& patchName : pressureBCs.sectionNames())
        {
            const auto& patchBC = pressureBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            Scalar value = S(0.0);
            if (bcType == "fixedValue")
            {
                value = patchBC.lookup<Scalar>("value");
            }
            else if (bcType == "fixedGradient")
            {
                value = patchBC.lookup<Scalar>("gradient");
            }

            registerBC(bcs, patchName, Field::p, bcType, value, "p");

            const auto& pType = *bcs[patchName][Field::p];

            hasFixedPressure = hasFixedPressure || pType.fixesValue();

            // Derive the p' boundary condition from p: fixed p becomes
            // p' = 0, zero-gradient p stays zero-gradient
            if (pType.fixesValue())
            {
                bcs[patchName][Field::pCorr] =
                    std::make_unique<FixedValue>(S(0.0));
            }
            else
            {
                bcs[patchName][Field::pCorr] =
                    std::make_unique<ZeroGradient>();
            }
        }
    }

    if (!hasFixedPressure)
    {
        Warning
        (
            "No fixedValue pressure boundary condition found. "
            "The pressure field has no reference value, which "
            "may cause a singular pressure matrix."
        );
    }

    // Inlet k per patch for omega's 'calculated' entries: a fixed value
    // carries itself, anything else falls back to the configured inlet k
    std::map<Name, Scalar> resolvedK;

    if (BCs.hasSection("k"))
    {
        const auto& kBCs = BCs.section("k");

        for (const auto& patchName : kBCs.sectionNames())
        {
            const auto& patchBC = kBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            // 'calculated' resolves loader-side: it needs the case config
            if (bcType == "fixedValue")
            {
                const Token valStr = patchBC.lookup<Token>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    value =
                        kOmegaSST::inletK
                        (
                            config.initialVelocity,
                            config.turbulenceIntensity
                        );

                    std::cout << std::format
                    (
                        "Inlet turbulence kinetic energy : {}\n",
                        value
                    );
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                registerBC(bcs, patchName, Field::k, bcType, value, "k");
                resolvedK[patchName] = value;
                continue;
            }

            Scalar value = S(0.0);
            if (bcType == "fixedGradient")
            {
                value = patchBC.lookup<Scalar>("gradient");
            }

            registerBC(bcs, patchName, Field::k, bcType, value, "k");
            resolvedK[patchName] =
                kOmegaSST::inletK
                (
                    config.initialVelocity,
                    config.turbulenceIntensity
                );
        }
    }

    if (BCs.hasSection("omega"))
    {
        const auto& omegaBCs = BCs.section("omega");

        for (const auto& patchName : omegaBCs.sectionNames())
        {
            const auto& patchBC = omegaBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            // 'calculated' resolves loader-side from the patch's inlet k
            if (bcType == "fixedValue")
            {
                const Token valStr = patchBC.lookup<Token>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    const auto kIterator = resolvedK.find(patchName);

                    if (kIterator == resolvedK.end())
                    {
                        FatalError
                        (
                            "Boundary condition not found for patch "
                          + patchName + " and field k"
                        );
                    }

                    value =
                        kOmegaSST::inletOmega
                        (
                            kIterator->second,
                            config.hydraulicDiameter
                        );

                    std::cout << std::format
                    (
                        "Inlet specific dissipation : {}\n",
                        value
                    );
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                registerBC(bcs, patchName, Field::omega, bcType, value, "omega");
                continue;
            }

            Scalar value = S(0.0);
            if (bcType == "fixedGradient")
            {
                value = patchBC.lookup<Scalar>("gradient");
            }

            registerBC(bcs, patchName, Field::omega, bcType, value, "omega");
        }
    }

    if (BCs.hasSection("nut"))
    {
        const auto& nutBCs = BCs.section("nut");

        for (const auto& patchName : nutBCs.sectionNames())
        {
            const auto& patchBC = nutBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            Scalar value = S(0.0);
            if (bcType == "fixedValue")
            {
                value = patchBC.lookup<Scalar>("value");
            }
            else if (bcType == "fixedGradient")
            {
                value = patchBC.lookup<Scalar>("gradient");
            }

            registerBC(bcs, patchName, Field::nut, bcType, value, "nut");
        }
    }

    BoundaryConditions bcManager(std::move(bcs), mesh);

    validateWallFunctionSetup(mesh, bcManager, config);

    if (config.debug)
    {
        bcManager.printSummary();
    }

    std::cout << std::format
    (
        "Boundary conditions set for {} patches.\n",
        mesh.patches().size()
    );

    return bcManager;
}

} // namespace BCLoader
