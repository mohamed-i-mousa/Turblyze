/******************************************************************************
 * @file BoundaryConditions.hpp
 * @brief Manages boundary conditions for the CFD solver
 *
 * @details This class provides functionality to set up, store, and apply
 * boundary conditions for different fields on mesh patches. It supports
 * various boundary condition types including fixed values, gradients, and
 * special conditions like no-slip.
 *
 * @class BoundaryConditions
 *
 * @example: usage hierarchy:
 * - Patch "inlet": U → fixed value, p → zero gradient
 * - Patch "outlet": U → zero gradient, p → fixed value
 * - Patch "wall": U → no-slip, p → zero gradient
 *****************************************************************************/

#pragma once

#include <vector>
#include <string>
#include <map>
#include <optional>

#include "Scalar.hpp"
#include "Vector.hpp"
#include "Face.hpp"
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
    void addPatch(BoundaryPatch patch);

    /**
     * @brief Set generic boundary condition
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field (U, p, k, omega, etc.)
     * @param bcData Boundary condition data
     * @return True if successfully set
     */
    bool setBC
    (
        const std::string& patchName,
        const std::string& fieldName,
        BoundaryData bcData
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
     * @brief Set wall function boundary condition type
     * @param patchName Name of the boundary patch
     * @param fieldName Name of the field
     * @param wallType The wall function type (K_WALL_FUNCTION, OMEGA_WALL_FUNCTION, or NUT_WALL_FUNCTION)
     * @return True if successfully set
     */
    bool setWallFunctionType
    (
        const std::string& patchName,
        const std::string& fieldName,
        BCType wallType
    );

// Accessor methods

    /**
     * @brief Get boundary patch by name
     * @param name Name of the patch to retrieve
     * @return Pointer to the found patch
     * @note Terminates the program if patch not found
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
     * @param fieldName Name of the field
     * @param phi Scalar field
     * @param face Boundary face
     * @param componentIdx Optional vector component index for extracting
     *        a scalar from a vector BC (0=x, 1=y, 2=z)
     * @return Boundary value based on boundary condition
     */
    [[nodiscard("Computed boundary face value needed for discretization")]]
    Scalar calculateBoundaryFaceValue
    (
        const std::string& fieldName,
        const ScalarField& phi,
        const Face& face,
        std::optional<int> componentIdx = std::nullopt
    ) const;

    /**
     * @brief Calculate boundary face value for vector field
     * @param fieldName Name of the field (e.g., "U")
     * @param phi Vector field
     * @param face Boundary face
     * @return Boundary vector value based on boundary condition
     * @note Defaults to zero-gradient if no BC is specified for the face
     */
    [[nodiscard("Computed boundary face value needed for discretization")]]
    Vector calculateBoundaryVectorFaceValue
    (
        const std::string& fieldName,
        const VectorField& phi,
        const Face& face
    ) const;

    /**
     * @brief Link boundary faces to their owning patches
     * @param faces All mesh faces (boundary faces get patch pointers set)
     * @note This is a state-changing initialization method; must be called
     *       exactly once before solving, and prevents further addPatch() calls
     */
    void linkFaces(std::vector<Face>& faces);

    /**
     * @brief Convert boundary condition type to string
     * @param bctype Boundary condition type enumeration
     * @return String representation of BC type
     * @note Terminates the program if unknown BC type
     * @note It's used in printSummary
     */
    static std::string bcTypeToString(BCType bctype);

    /**
     * @brief Validate that all patch names in boundary data
     *        exist in the mesh patches
     * @note Terminates the program if an unknown patch name is found
     */
    void validatePatchNames() const;

    /// Print summary of all boundary conditions
    void printSummary() const;

private:

// Private members

    /// Nested map: patch name → field name → boundary data
    std::map<std::string, std::map<std::string, BoundaryData>>
    patchBoundaryData_;

    /// Vector of all boundary patches
    std::vector<BoundaryPatch> patches_;

    /// True after linkFaces() — prevents addPatch() after linking
    bool linked_ = false;
};
