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

#include <iostream>
#include <string>
#include <iomanip>
#include <chrono>

#include "Scalar.h"
#include "CFDApplication.h"


int main(int argc, char* argv[])
{
    // Start timing the total execution
    auto startTime = std::chrono::high_resolution_clock::now();

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
)" << std::endl;

    std::cout
        << "Running with precision: " << SCALAR_MODE
        << std::endl;

    std::cout 
        << std::fixed << std::setprecision(6);

    std::string caseFile = "../defaultCase";

    if (argc > 1)
    {
        caseFile = argv[1];

        std::cout
            << "Using case file: " << caseFile
            << std::endl;
    }
    else
    {
        std::cout
            << "Using default case: " << caseFile
            << std::endl;
    }

    CFDApplication app(caseFile);
    app.run();

    auto endTime = std::chrono::high_resolution_clock::now();

    auto duration = 
        std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);

    std::cout
        << "\n--- Simulation Complete ---" << std::endl;

    std::cout
        << "Total execution time: " << duration.count()
        << " seconds" << std::endl;

    return 0;
}