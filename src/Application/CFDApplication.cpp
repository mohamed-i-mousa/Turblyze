/******************************************************************************
 * @file CFDApplication.cpp
 * @brief Top-level application orchestrator for the CFD solver
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "CFDApplication.h"

// Standard library headers
#include <iostream>

// External library headers
#include <eigen3/Eigen/Core>
#include <omp.h>

// Project headers
#include "BoundaryConditionsLoader.h"
#include "BoundaryConditions.h"
#include "CaseConfiguration.h"
#include "CaseReader.h"
#include "MeshCreator.h"
#include "PostProcess.h"
#include "SolverSetup.h"
#include "SIMPLE.h"

// *************************** Internal helpers *******************************

namespace
{

void initParallelism(int numThreads)
{
    omp_set_num_threads(numThreads);
    Eigen::setNbThreads(numThreads);

    std::cout
        << "OpenMP threads: " << numThreads << '\n';
}

} // namespace

// ************************* Special Member Functions *************************

CFDApplication::CFDApplication(const std::string& caseFilePath)
:
    caseFilePath_(caseFilePath)
{}

CFDApplication::~CFDApplication() noexcept = default;

// *********************************** run ***********************************

void CFDApplication::run()
{
    std::cout
        << '\n' << "--- 0. Loading Case ---" << '\n';

    // Read the case file and load configuration
    CaseReader caseReader(caseFilePath_);
    const CaseConfiguration config = CaseConfig::loadConfiguration(caseReader);

    // Initialize parallelism
    initParallelism(config.numThreads);

    // Create mesh
    Mesh mesh = MeshCreator::create(config);

    // Load boundary conditions
    BoundaryConditions bcManager;
    BCLoader::load(caseReader, config, mesh, bcManager);

    // Configure solver
    SolverModules modules;
    SolverSetup::configure(modules, mesh, bcManager, config);
    SolverSetup::logSetup(modules, config);

    std::cout
        << '\n' << "--- 4. Solving Steady-State Flow with SIMPLE ---"
        << '\n';

    // Run the solver
    modules.solver->solve();

    // Post-process results
    PostProcess::reportStatistics(*modules.solver, config);
    PostProcess::exportResults
    (
        *modules.solver,
        modules.turbulenceModel.get(),
        mesh,
        config
    );
}
