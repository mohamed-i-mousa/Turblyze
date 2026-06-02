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

// ********************************** Headers *********************************

// Standard library headers
#include <vector>
#include <string>
#include <map>
#include <optional>
#include <span>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "BoundaryPatch.h"
#include "BoundaryData.h"
#include "CellData.h"
#include "Field.h"

// ************************* class BoundaryConditions *************************

class BoundaryConditions
{
public:

// ****************************** Setter Methods ******************************

    /// Add a boundary patch from mesh reader
    void addPatch(BoundaryPatch patch);

    /// Set generic boundary condition
    void setBC
    (
        const std::string& patchName,
        Field field,
        BoundaryData bcData
    )
    {
        patchBoundaryData_[patchName][field] = std::move(bcData);
    }

    /// Set fixed scalar value boundary condition
    void setFixedValue
    (
        const std::string& patchName,
        Field field,
        Scalar value
    );

    /// Set fixed scalar gradient boundary condition
    void setFixedGradient
    (
        const std::string& patchName,
        Field field,
        Scalar gradient
    );

    /// Set zero gradient boundary condition
    void setZeroGradient
    (
        const std::string& patchName,
        Field field
    );

    /// Set no-slip boundary condition
    void setNoSlip
    (
        const std::string& patchName,
        Field field
    );

    /// Set wall function boundary condition type
    void setWallFunctionType
    (
        const std::string& patchName,
        Field field,
        BCType wallType
    );

// ***************************** Accessor Methods *****************************

    /// Get boundary patch by name
    [[nodiscard]] const BoundaryPatch& patch(const std::string& name) const;

    /// Get all boundary patches
    [[nodiscard]] const std::vector<BoundaryPatch>& patches() const noexcept
    {
        return patches_;
    }

    /// Get number of patches
    [[nodiscard]] size_t numPatches() const noexcept
    {
        return patches_.size();
    }

    /// Get boundary condition for a field on a patch
    [[nodiscard]] const BoundaryData& fieldBC
    (
        const std::string& patchName,
        Field field
    ) const;

    /// Calculate boundary face value for scalar field
    [[nodiscard]] Scalar boundaryFaceValue
    (
        Field field,
        const ScalarField& phi,
        const Face& face
    ) const;

    /// Link boundary faces to their owning patches
    void linkFaces(std::span<Face> faces);

    /// Validate boundary condition patch names against mesh patch names
    void validatePatchNames() const;

    /// Print summary of all boundary conditions
    void printSummary() const;

// ****************************** Private Members *****************************

private:

    /// Nested map: patch name → field → boundary data
    std::map<std::string, std::map<Field, BoundaryData>>
    patchBoundaryData_;

    /// Vector of all boundary patches
    std::vector<BoundaryPatch> patches_;

    /// True after linkFaces() — prevents addPatch() after linking
    bool linked_ = false;
};
