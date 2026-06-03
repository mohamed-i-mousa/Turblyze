/******************************************************************************
 * @file main.cpp
 * @brief Main entry point for the 3D incompressible CFD solver
 *
 * This file contains the main function that launches the CFD simulation
 * via the CFDApplication driver class. It handles command-line argument
 * parsing, and timing
 *
 * @author Mohamed Mousa
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <iostream>
#include <iomanip>
#include <chrono>

// Project headers
#include "Scalar.h"
#include "StringTypes.h"
#include "CFDApplication.h"

// *********************************** main ***********************************

int main(int argc, char* argv[])
{
    // Start timing the total execution
    const auto startTime = std::chrono::high_resolution_clock::now();

    std::cout << R"(
  ~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~··~·~

  ████████╗██╗   ██╗██████╗ ██████╗ ██╗  ██╗   ██╗███████╗███████╗
  ╚══██╔══╝██║   ██║██╔══██╗██╔══██╗██║  ╚██╗ ██╔╝╚══███╔╝██╔════╝
     ██║   ██║   ██║██████╔╝██████╔╝██║   ╚████╔╝   ███╔╝ █████╗
     ██║   ██║   ██║██╔══██╗██╔══██╗██║    ╚██╔╝   ███╔╝  ██╔══╝
     ██║   ╚██████╔╝██║  ██║██████╔╝███████╗██║   ███████╗███████╗
     ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═════╝ ╚══════╝╚═╝   ╚══════╝╚══════╝

           3D Incompressible Navier-Stokes Solver v1.0

  ~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~··~·~
)" << '\n';

    std::cout
        << "Running with precision: " << SCALAR_MODE << '\n';

    std::cout 
        << std::fixed << std::setprecision(6);

    FilePath caseFile = "../defaultCase";

    if (argc > 1)
    {
        caseFile = argv[1];

        std::cout
            << "Using case file: " << caseFile << '\n';
    }
    else
    {
        std::cout
            << "Using default case: " << caseFile << '\n';
    }

    CFDApplication app(caseFile);
    app.run();

    const auto endTime = std::chrono::high_resolution_clock::now();

    const auto duration =
        std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);

    std::cout
        << '\n' << "--- Simulation Complete ---" << '\n';

    std::cout
        << "Total execution time: " << duration.count() << " seconds" << '\n'
        << '\n';

    return 0;
}
