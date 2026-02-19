/******************************************************************************
 * @file BoundaryConditions.hpp
 * @brief Manages boundary conditions for the CFD solver
 * 
 * This class provides functionality to set up, store, and apply boundary
 * conditions for different fields on mesh patches. It supports various
 * boundary condition types including fixed values, gradients, and
 * special conditions like no-slip.
 *
 * @class BoundaryConditions
 *
 * @example: usage hierarchy:
 * - Patch "inlet": U → fixed value, p → zero gradient
 * - Patch "outlet": U → zero gradient, p → fixed value  
 * - Patch "wall": U → no-slip, p → zero gradient
 *****************************************************************************/

#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <vector>
#include <string>
#include <map>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "BoundaryPatch.hpp"
#include "BoundaryData.hpp"
#include "CellData.hpp"

class BoundaryConditions 
{
public:

    /// Default constructor
    BoundaryConditions() = default;

// Setter methods

    /**
     * @brief Add a boundary patch from mesh reader
     * @param patch Boundary patch to add
     */
    void addPatch(const BoundaryPatch& patch);

    /**
     * @brief Set generic boundary condition
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field (U, p, k, omega, etc.)
     * @param BCSetup Boundary condition setup
     * @return True if successfully set
     */
    bool setBC
    (
        const std::string& patchName,
        const std::string& fieldName,
        BoundaryData BCSetup
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

    /// Set OpenFOAM-like kqR wall function boundary condition
    bool setKWallFunction
    (
        const std::string& patchName,
        const std::string& fieldName
    );

    /// Set OpenFOAM-like omega wall function boundary condition
    bool setOmegaWallFunction
    (
        const std::string& patchName,
        const std::string& fieldName
    );

    /// Set OpenFOAM-like nutk wall function boundary condition
    bool setNutWallFunction
    (
        const std::string& patchName,
        const std::string& fieldName
    );

// Accessor methods

    /**
     * @brief Get boundary patch by name
     * @param name Name of the patch to retrieve
     * @return Pointer to patch, or nullptr if not found
     * @throws std::runtime_error if patch not found
     */
    const BoundaryPatch* patch(const std::string& name) const;

    /** 
     * @brief Get all boundary patches 
     * @return Const reference to vector of boundary patches 
     */
    const std::vector<BoundaryPatch>& patches() const { return patches_; }

    /** 
     * @brief Get number of patches 
     * @return Number of boundary patches 
     */
    size_t numPatches() const { return patches_.size(); }

    /**
     * @brief Get boundary condition for a field on a patch
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field
     * @return Pointer to boundary data, or nullptr if not found
     */
    const BoundaryData* fieldBC
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
     * @param fieldName Name of the field (e.g., "U")
     * @return Boundary vector value based on boundary condition
     * @throws std::runtime_error if face not found in boundary patches
     */
    Vector calculateBoundaryVectorFaceValue
    (
        const Face& face,
        const VectorField& phi,
        const std::string& fieldName
    ) const;

    /**
     * @brief Link boundary faces to their owning patches
     * @param faces All mesh faces (boundary faces get patch pointers set)
     */
    void linkFaces(std::vector<Face>& faces) const;

    /**
     * @brief Convert boundary condition type to string
     * @param bctype Boundary condition type enumeration
     * @return String representation of BC type
     * @throws std::runtime_error if unknown BC type
     * @note It's used in printSummary
     */
    std::string bcTypeToString(BCType bctype) const;

    /**
     * @brief Print summary of all boundary conditions
     */
    void printSummary() const;

private:

// Private members

    /// Nested map: patch name → field name → boundary data
    std::map<std::string, std::map<std::string, BoundaryData>>
    patchBoundaryData_;
    
    /// Vector of all boundary patches
    std::vector<BoundaryPatch> patches_;
};

#endif // BOUNDARY_CONDITIONS_HPP
