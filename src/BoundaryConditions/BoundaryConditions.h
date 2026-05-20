/******************************************************************************
 * @file BoundaryConditions.h
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
#include <span>

#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "BoundaryPatch.h"
#include "BoundaryData.h"
#include "CellData.h"
#include "Field.h"


class BoundaryConditions
{
public:

// Setter methods

    /**
     * @brief Add a boundary patch from mesh reader
     * @param patch Boundary patch to add
     */
    void addPatch(BoundaryPatch patch);

    /**
     * @brief Set generic boundary condition
     * @param patchName Name of the boundary patch
     * @param field Field identifier
     * @param bcData Boundary condition data
     */
    void setBC
    (
        const std::string& patchName,
        Field field,
        BoundaryData bcData
    )
    {
        patchBoundaryData_[patchName][field] = std::move(bcData);
    }

    /**
     * @brief Set fixed scalar value boundary condition
     * @param patchName Name of the boundary patch
     * @param field Field identifier
     * @param value Scalar value to fix at boundary
     */
    void setFixedValue
    (
        const std::string& patchName,
        Field field,
        Scalar value
    );

    /**
     * @brief Set fixed scalar gradient boundary condition
     * @param patchName Name of the boundary patch
     * @param field Field identifier
     * @param gradient Scalar gradient to fix at boundary
     */
    void setFixedGradient
    (
        const std::string& patchName,
        Field field,
        Scalar gradient
    );

    /**
     * @brief Set zero gradient boundary condition
     * @param patchName Name of the boundary patch
     * @param field Field identifier
     */
    void setZeroGradient
    (
        const std::string& patchName,
        Field field
    );

    /**
     * @brief Set no-slip boundary condition (for velocity)
     * @param patchName Name of the boundary patch
     * @param field Field identifier
     */
    void setNoSlip
    (
        const std::string& patchName,
        Field field
    );

    /**
     * @brief Set wall function boundary condition type
     * @param patchName Name of the boundary patch
     * @param field Field identifier
     * @param wallType The wall function type (K_WALL_FUNCTION,
     *                 OMEGA_WALL_FUNCTION, or NUT_WALL_FUNCTION)
     */
    void setWallFunctionType
    (
        const std::string& patchName,
        Field field,
        BCType wallType
    );

// Accessor methods

    /**
     * @brief Get boundary patch by name
     * @param name Name of the patch to retrieve
     * @return Reference to the found patch
     * @note Terminates the program if patch not found
     */
    [[nodiscard]] const BoundaryPatch& patch(const std::string& name) const;

    /**
     * @brief Get all boundary patches
     * @return Const reference to vector of boundary patches
     */
    [[nodiscard]] const std::vector<BoundaryPatch>& patches() const noexcept
    {
        return patches_;
    }

    /**
     * @brief Get number of patches
     * @return Number of boundary patches
     */
    [[nodiscard]] size_t numPatches() const noexcept
    {
        return patches_.size();
    }

    /**
     * @brief Get boundary condition for a field on a patch
     * @param patchName Name of the boundary patch
     * @param field Field identifier
     * @return Reference to boundary data
     * @note Terminates the program if not found
     */
    [[nodiscard]] const BoundaryData& fieldBC
    (
        const std::string& patchName,
        Field field
    ) const;

    /**
     * @brief Calculate boundary face value for scalar field
     * @param field Field identifier
     * @param phi Scalar field
     * @param face Boundary face
     * @return Boundary value based on boundary condition
     */
    [[nodiscard]] Scalar boundaryFaceValue
    (
        Field field,
        const ScalarField& phi,
        const Face& face
    ) const;

    /**
     * @brief Link boundary faces to their owning patches
     * @param faces All mesh faces (boundary faces get patch pointers set)
     * @note This is a state-changing initialization method; must be called
     *       exactly once before solving, and prevents further addPatch() calls
     */
    void linkFaces(std::span<Face> faces);

    /**
     * @brief Convert boundary condition type to string
     * @param bctype Boundary condition type enumeration
     * @return String representation of BC type
     * @note Terminates the program if unknown BC type
     * @note It's used in printSummary
     */
    [[nodiscard]] static std::string bcTypeToString(BCType bctype) noexcept;

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

    /// Nested map: patch name → field → boundary data
    std::map<std::string, std::map<Field, BoundaryData>>
    patchBoundaryData_;

    /// Vector of all boundary patches
    std::vector<BoundaryPatch> patches_;

    /// True after linkFaces() — prevents addPatch() after linking
    bool linked_ = false;
};
