/******************************************************************************
 * @file PostProcess.h
 * @brief After-solve reporting and result export
 *****************************************************************************/

#pragma once

// *************************** Forward Declarations ***************************

class kOmegaSST;
class Mesh;
class SIMPLE;
struct CaseConfiguration;

// *************************** namespace PostProcess **************************

namespace PostProcess
{

/// Extract solution fields and print flow statistics
void reportStatistics
(
    const SIMPLE& solver,
    const CaseConfiguration& config
);

/// Write VTK volume and wall-boundary results
void exportResults
(
    const SIMPLE& solver,
    const kOmegaSST* turbulence,
    const Mesh& mesh,
    const CaseConfiguration& config
);

} // namespace PostProcess
