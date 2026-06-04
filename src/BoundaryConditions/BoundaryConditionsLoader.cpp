/******************************************************************************
 * @file BoundaryConditionsLoader.cpp
 * @brief Case-file boundary condition registration
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "BoundaryConditionsLoader.h"

// Standard library headers
#include <iostream>

// Project headers
#include "CaseReader.h"
#include "ErrorHandler.h"
#include "Logger.h"
#include "kOmegaSST.h"

// **************************** namespace BCLoader ****************************

namespace BCLoader
{

// ***************************** Internal Helpers *****************************

namespace
{

[[noreturn]] void unknownBCType
(
    NameRef bcType,
    NameRef fieldName,
    NameRef patchName,
    MessageRef validList
)
{
    FatalError
    (
        "Unknown boundary condition type '" + Name(bcType)
      + "' for field '" + Name(fieldName)
      + "' on patch '" + Name(patchName)
      + "'. Valid types: " + Message(validList)
    );
}


void validateWallFunctionSetup
(
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    if (!config.turbulenceEnabled)
    {
        return;
    }

    for (const auto& patch : mesh.patches())
    {
        if (patch.type() != PatchType::wall)
        {
            continue;
        }

        const Name& patchName = patch.patchName();

        const bool kIsWF =
            bcManager.fieldBC(patchName, Field::k).type()
         == BCType::kWallFunction;
        const bool omegaIsWF =
            bcManager.fieldBC(patchName, Field::omega).type()
         == BCType::omegaWallFunction;
        const bool nutIsWF =
            bcManager.fieldBC(patchName, Field::nut).type()
         == BCType::nutWallFunction;

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

} // namespace (unnamed)


// *********************************** Load ***********************************

void load
(
    const CaseReader& reader,
    const CaseConfiguration& config,
    Mesh& mesh,
    BoundaryConditions& bcManager
)
{
    std::cout << '\n';
    Logger::sectionHeader("Setting Boundary Conditions");

    for (const auto& patch : mesh.patches())
    {
        bcManager.addPatch(patch);
    }

    bcManager.linkFaces(mesh.faces());

    for (const auto& face : mesh.faces())
    {
        if (face.isBoundary() && !face.patch().has_value())
        {
            FatalError
            (
                "Boundary face "
              + std::to_string(face.idx())
              + " has no patch after linking."
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

            if (bcType == "fixedValue")
            {
                const Vector value = patchBC.lookup<Vector>("value");
                bcManager.setFixedValue(patchName, Field::Ux, value.x());
                bcManager.setFixedValue(patchName, Field::Uy, value.y());
                bcManager.setFixedValue(patchName, Field::Uz, value.z());
            }
            else if (bcType == "noSlip")
            {
                bcManager.setNoSlip(patchName, Field::Ux);
                bcManager.setNoSlip(patchName, Field::Uy);
                bcManager.setNoSlip(patchName, Field::Uz);
            }
            else if (bcType == "zeroGradient")
            {
                bcManager.setZeroGradient(patchName, Field::Ux);
                bcManager.setZeroGradient(patchName, Field::Uy);
                bcManager.setZeroGradient(patchName, Field::Uz);
            }
            else
            {
                unknownBCType
                (
                    bcType, "U", patchName,
                    "fixedValue, noSlip, zeroGradient"
                );
            }
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

            if (bcType == "fixedValue")
            {
                const Scalar value = patchBC.lookup<Scalar>("value");
                bcManager.setFixedValue(patchName, Field::p, value);
                hasFixedPressure = true;
            }
            else if (bcType == "zeroGradient")
            {
                bcManager.setZeroGradient(patchName, Field::p);
            }
            else
            {
                unknownBCType
                (
                    bcType, "p", patchName,
                    "fixedValue, zeroGradient"
                );
            }
        }

        for (const auto& patchName : pressureBCs.sectionNames())
        {
            const auto& patchBC = pressureBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            if (bcType == "fixedValue")
            {
                bcManager.setFixedValue(patchName, Field::pCorr, S(0.0));
            }
            else if (bcType == "zeroGradient")
            {
                bcManager.setZeroGradient(patchName, Field::pCorr);
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

    if (BCs.hasSection("k"))
    {
        const auto& kBCs = BCs.section("k");

        for (const auto& patchName : kBCs.sectionNames())
        {
            const auto& patchBC = kBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            if (bcType == "fixedValue")
            {
                const Token valStr =
                    patchBC.lookup<Token>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    value =
                        kOmegaSST::inletK
                        (
                            config.initialVelocity,
                            config.turbulenceIntensity
                        );

                    std::cout
                        << "Inlet turbulence kinetic energy : " << value
                        << '\n';
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                bcManager.setFixedValue(patchName, Field::k, value);
            }
            else if (bcType == "kWallFunction")
            {
                bcManager.setWallFunctionType
                (
                    patchName,
                    Field::k,
                    BCType::kWallFunction
                );
            }
            else if (bcType == "zeroGradient")
            {
                bcManager.setZeroGradient(patchName, Field::k);
            }
            else
            {
                unknownBCType
                (
                    bcType, "k", patchName,
                    "fixedValue, kWallFunction, zeroGradient"
                );
            }
        }
    }

    if (BCs.hasSection("omega"))
    {
        const auto& omegaBCs = BCs.section("omega");

        for (const auto& patchName : omegaBCs.sectionNames())
        {
            const auto& patchBC = omegaBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            if (bcType == "fixedValue")
            {
                const Token valStr =
                    patchBC.lookup<Token>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    const BoundaryData& kPatchBC =
                        bcManager.fieldBC(patchName, Field::k);

                    const Scalar kValue =
                        kPatchBC.type() == BCType::fixedValue
                      ? kPatchBC.fixedScalarValue()
                      : kOmegaSST::inletK
                        (
                            config.initialVelocity,
                            config.turbulenceIntensity
                        );

                    value =
                        kOmegaSST::inletOmega
                        (
                            kValue,
                            config.hydraulicDiameter
                        );

                    std::cout
                        << "Inlet specific dissipation : " << value << '\n';
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                bcManager.setFixedValue(patchName, Field::omega, value);
            }
            else if (bcType == "omegaWallFunction")
            {
                bcManager.setWallFunctionType
                (
                    patchName,
                    Field::omega,
                    BCType::omegaWallFunction
                );
            }
            else if (bcType == "zeroGradient")
            {
                bcManager.setZeroGradient(patchName, Field::omega);
            }
            else
            {
                unknownBCType
                (
                    bcType, "omega", patchName,
                    "fixedValue, omegaWallFunction, zeroGradient"
                );
            }
        }
    }

    if (BCs.hasSection("nut"))
    {
        const auto& nutBCs = BCs.section("nut");

        for (const auto& patchName : nutBCs.sectionNames())
        {
            const auto& patchBC = nutBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            if (bcType == "fixedValue")
            {
                const Scalar value = patchBC.lookup<Scalar>("value");
                bcManager.setFixedValue(patchName, Field::nut, value);
            }
            else if (bcType == "zeroGradient")
            {
                bcManager.setZeroGradient(patchName, Field::nut);
            }
            else if (bcType == "nutWallFunction")
            {
                bcManager.setWallFunctionType
                (
                    patchName,
                    Field::nut,
                    BCType::nutWallFunction
                );
            }
            else
            {
                unknownBCType
                (
                    bcType, "nut", patchName,
                    "fixedValue, zeroGradient, nutWallFunction"
                );
            }
        }
    }

    bcManager.validatePatchNames();

    validateWallFunctionSetup(mesh, bcManager, config);

    if (config.debug)
    {
        bcManager.printSummary();
    }

    std::cout
        << "Boundary conditions set for "
        << mesh.patches().size() << " patches." << '\n';
}

} // namespace BCLoader
