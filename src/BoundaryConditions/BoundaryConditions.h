#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <vector>
#include <string>
#include <map>

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"
#include "BoundaryData.h"
#include "CellData.h"

class BoundaryConditions {
public:
    // ----- Constructor ----- //
    BoundaryConditions() = default;

    // ----- Methods ----- //

    // Add patch from the mesh reader
    void addPatch(const BoundaryPatch& patch);

    // Get a patch by its name
    const BoundaryPatch* getPatch(const std::string& name) const;

    // ----- Boundary Condition Setters ----- //
    // Generic setter
    bool setBC(const std::string& patchName, const std::string& fieldName, BoundaryData bc_config);

    // Fixed value setters
    bool setFixedValue(const std::string& patchName, const std::string& fieldName, Scalar value);
    bool setFixedValue(const std::string& patchName, const std::string& fieldName, const Vector& value);

    // Gradient setters
    bool setFixedGradient(const std::string& patchName, const std::string& fieldName, Scalar gradient);
    bool setFixedGradient(const std::string& patchName, const std::string& fieldName, const Vector& gradient);

    // Special boundary condition setters
    bool setZeroGradient(const std::string& patchName, const std::string& fieldName);
    bool setNoSlip(const std::string& patchName, const std::string& fieldName);
    bool setSymmetry(const std::string& patchName, const std::string& fieldName);

    // ----- Boundary Condition Getters ----- //
    // Get boundary condition configuration
    const BoundaryData* getFieldBC(const std::string& patchName, const std::string& fieldName) const;

    // ----- Boundary Face Value Calculations ----- //
    // Calculate boundary face value based on boundary conditions
    Scalar calculateBoundaryFaceValue(
        const Face& face,
        const ScalarField& phi,
        const std::string& fieldName) const;

    // Calculate boundary face vector value based on boundary conditions
    Vector calculateBoundaryFaceVectorValue(
        const Face& face,
        const VectorField& phi,
        const std::string& fieldName) const;

    // ----- Utility & Debugging Methods ----- //
    // Convert BCType enum to string for printing
    std::string bcTypeToString(BCType bctype) const;

    // Print boundary conditions data 
    void printSummary() const;

private:
    // ----- Private Member Variables ----- //
    // Shared cache for face-to-patch mapping
    mutable std::map<size_t, const BoundaryPatch*> faceToPatchCache;
    mutable bool cacheBuilt = false;
    
    // ----- Private Helper Methods ----- //
    // Helper method to ensure cache is built
    void ensureFaceToPatchCacheBuilt() const;

public:
    // ----- Public Member Variables ----- //
    // The main storage of boundary conditions setup
    /* Example:
    bcManager
    ├── "inlet"
    │   ├── "U"     --> BoundaryData(FIXED_VALUE, VECTOR, (10,0,0))
    │   ├── "p"     --> BoundaryData(ZERO_GRADIENT)
    │   ├── "k"     --> BoundaryData(FIXED_VALUE, SCALAR, 0.1)
    │   └── "omega" --> BoundaryData(FIXED_VALUE, SCALAR, 1.0)
    │
    ├── "outlet"
    │   ├── "U"     --> BoundaryData(ZERO_GRADIENT)
    │   ├── "p"     --> BoundaryData(FIXED_VALUE, SCALAR, 101325.0)
    │   ├── "k"     --> BoundaryData(ZERO_GRADIENT)
    │   └── "omega" --> BoundaryData(ZERO_GRADIENT)
    */
    std::map<std::string, std::map<std::string, BoundaryData>> patchBoundaryData;
    std::vector<BoundaryPatch> patches;
};

#endif