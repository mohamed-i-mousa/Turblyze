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

#include "Scalar.hpp"
#include "CFDApplication.hpp"


int main(int argc, char* argv[])
{
    // Start timing the total execution
    auto startTime = std::chrono::high_resolution_clock::now();

    std::cout
        << "--- Welcome to the 3D Incompressible CFD Solver ---"
        << std::endl;

    std::cout
        << "Running with precision: " << SCALAR_MODE
        << std::endl;

    std::cout 
        << std::fixed << std::setprecision(6);

    try
    {
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
    }
    catch (const std::exception& e)
    {
        std::cerr
            << std::endl << "Error: " << e.what() << std::endl;
        return 1;
    }

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