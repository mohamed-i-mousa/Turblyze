/******************************************************************************
 * @file CFDApplication.h
 * @brief Top-level application orchestrator for the CFD solver
 *
 * @details CFDApplication owns only the case-file path and coordinates the
 * simulation phases in run(). Phase-specific parsing, setup, solver assembly,
 * reporting, and export live in focused modules.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <string>

// *************************** class CFDApplication ***************************

class CFDApplication
{
public:

// ************************* Special Member Functions *************************

    /// Constructor for CFDApplication
    explicit CFDApplication(const std::string& caseFilePath);

    /// Copy constructor and assignment - Not copyable
    CFDApplication(const CFDApplication&) = delete;
    CFDApplication& operator=(const CFDApplication&) = delete;

    /// Move constructor and assignment - Not movable
    CFDApplication(CFDApplication&&) = delete;
    CFDApplication& operator=(CFDApplication&&) = delete;

    /// Destructor
    ~CFDApplication() noexcept;

// ******************************** Solver Run ********************************

    /// Run the full simulation
    void run();

// ****************************** Private Members *****************************

private:

    /// Path to case file
    std::string caseFilePath_;
};
