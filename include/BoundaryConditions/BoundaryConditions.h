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

/**
 * @brief Manages boundary conditions for the CFD solver
 * 
 * This class provides functionality to set up, store, and apply boundary
 * conditions for different fields on mesh patches. It supports various
 * boundary condition types including fixed values, gradients, and
 * special conditions like no-slip.
 * 
 * Example usage hierarchy:
 * - Patch "inlet": U → fixed value, p → zero gradient
 * - Patch "outlet": U → zero gradient, p → fixed value  
 * - Patch "wall": U → no-slip, p → zero gradient
 */
class BoundaryConditions 
{
public:
    /**
     * @brief Default constructor
     */
    BoundaryConditions() = default;

    /**
     * @brief Add a boundary patch from mesh reader
     * @param patch Boundary patch to add
     */
    void addPatch(const BoundaryPatch& patch);

    /**
     * @brief Get boundary patch by name
     * @param name Name of the patch to retrieve
     * @return Pointer to patch, or nullptr if not found
     * @throws std::runtime_error if patch not found
     */
    const BoundaryPatch* getPatch(const std::string& name) const;

    /**
     * @brief Set generic boundary condition
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field (U, p, k, omega, etc.)
     * @param bc_config Boundary condition configuration
     * @return True if successfully set
     */
    bool setBC
    (
        const std::string& patchName, 
        const std::string& fieldName, 
        BoundaryData bc_config
    );

    /**
     * @brief Set fixed scalar value boundary condition
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field
     * @param value Scalar value to fix at boundary
     * @return True if successfully set
     */
    bool setFixedValue
    (
        const std::string& patchName, 
        const std::string& fieldName, 
        Scalar value
    );
    
    /**
     * @brief Set fixed vector value boundary condition
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field
     * @param value Vector value to fix at boundary
     * @return True if successfully set
     */
    bool setFixedValue
    (
        const std::string& patchName, 
        const std::string& fieldName, 
        const Vector& value
    );

    /**
     * @brief Set fixed scalar gradient boundary condition
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field
     * @param gradient Scalar gradient to fix at boundary
     * @return True if successfully set
     */
    bool setFixedGradient
    (
        const std::string& patchName, 
        const std::string& fieldName, 
        Scalar gradient
    );
    
    /**
     * @brief Set fixed vector gradient boundary condition
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field
     * @param gradient Vector gradient to fix at boundary
     * @return True if successfully set
     */
    bool setFixedGradient
    (
        const std::string& patchName, 
        const std::string& fieldName, 
        const Vector& gradient
    );

    /**
     * @brief Set zero gradient boundary condition
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field
     * @return True if successfully set
     */
    bool setZeroGradient
    (
        const std::string& patchName, 
        const std::string& fieldName
    );
    
    /**
     * @brief Set no-slip boundary condition (for velocity)
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field
     * @return True if successfully set
     */
    bool setNoSlip
    (
        const std::string& patchName, 
        const std::string& fieldName
    );
    

    /**
     * @brief Get boundary condition for a field on a patch
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field
     * @return Pointer to boundary data, or nullptr if not found
     */
    const BoundaryData* getFieldBC
    (
        const std::string& patchName, 
        const std::string& fieldName
    ) const;

    /**
     * @brief Calculate boundary face value for scalar field
     * @param face Boundary face
     * @param phi Scalar field
     * @param fieldName Name of the field
     * @return Boundary value based on boundary condition
     * @throws std::runtime_error if face not found in boundary patches
     */
    Scalar calculateBoundaryFaceValue
    (
        const Face& face,
        const ScalarField& phi,
        const std::string& fieldName
    ) const;

    /**
     * @brief Calculate boundary face value for vector field
     * @param face Boundary face
     * @param phi Vector field
     * @param fieldName Name of the field
     * @return Boundary value based on boundary condition
     * @throws std::runtime_error if face not found in boundary patches
     */
    Vector calculateBoundaryFaceVectorValue
    (
        const Face& face,
        const VectorField& phi,
        const std::string& fieldName
    ) const;

    /**
     * @brief Convert boundary condition type to string
     * @param bctype Boundary condition type enumeration
     * @return String representation of BC type
     * @throws std::runtime_error if unknown BC type
     */
    std::string bcTypeToString(BCType bctype) const;

    /**
     * @brief Print summary of all boundary conditions
     */
    void printSummary() const;

    /**
     * @brief Get all boundary patches
     * @return Const reference to vector of boundary patches
     */
    const std::vector<BoundaryPatch>& getPatches() const { return patches; }

    /**
     * @brief Get number of patches
     * @return Number of boundary patches
     */
    size_t getNumPatches() const { return patches.size(); }

private:
    /// Nested map: patch name → field name → boundary data
    std::map<std::string, std::map<std::string, BoundaryData>> 
        patchBoundaryData;
    
    /// Vector of all boundary patches
    std::vector<BoundaryPatch> patches;
    /// Cache for fast face-to-patch mapping
    mutable std::map<size_t, const BoundaryPatch*> faceToPatchCache;
    
    /// Flag indicating if cache has been built
    mutable bool cacheBuilt = false;
    
    /**
     * @brief Ensure face-to-patch cache is built for efficient lookup
     */
    void ensureFaceToPatchCacheBuilt() const;
};

#endif