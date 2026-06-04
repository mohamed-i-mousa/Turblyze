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
#include "Logger.h"
#include "MeshCreator.h"
#include "PostProcess.h"
#include "SolverSetup.h"
#include "SIMPLE.h"

// ***************************** Internal Helpers *****************************

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

CFDApplication::CFDApplication(const FilePath& caseFile)
:
    caseFile_(caseFile)
{}

CFDApplication::~CFDApplication() noexcept = default;

// ******************************** Solver Run ********************************

void CFDApplication::run()
{
    std::cout << '\n';
    Logger::sectionHeader("Loading Case");

    // Read the case file and load configuration
    CaseReader caseReader(caseFile_);
    const CaseConfiguration config = CaseConfig::loadConfiguration(caseReader);

    // Initialize parallelism
    initParallelism(static_cast<int>(config.numThreads));

    // Create mesh
    Mesh mesh = MeshCreator::create(config);

    // Load boundary conditions
    BoundaryConditions bcManager;
    BCLoader::load(caseReader, config, mesh, bcManager);

    // Configure solver
    SolverModules modules;
    SolverSetup::configure(modules, mesh, bcManager, config);
    SolverSetup::logSetup(modules, config);

    // Run the solver (SIMPLE::solve prints its own framed banner)
    modules.solver->solve();

    // Post-process results
    PostProcess::reportStatistics(*modules.solver);
    PostProcess::exportResults
    (
        *modules.solver,
        *modules.turbulenceModel,
        mesh,
        config
    );
}
