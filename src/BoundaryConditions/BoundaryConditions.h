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
#include <map>
#include <optional>

// Project headers
#include "Scalar.h"
#include "MeshContainers.h"
#include "Face.h"
#include "BoundaryPatch.h"
#include "BoundaryData.h"
#include "CellData.h"
#include "Field.h"
#include "Integer.h"
#include "StringTypes.h"

// ************************* class BoundaryConditions *************************

class BoundaryConditions
{
public:

    using PatchBoundaryDataMap =
        std::map<Name, std::map<Field, BoundaryData>>;

// ****************************** Setter Methods ******************************

    /// Add a boundary patch from mesh reader
    void addPatch(BoundaryPatch patch);

    /// Set generic boundary condition
    void setBC
    (
        const Name& patchName,
        Field field,
        BoundaryData bcData
    )
    {
        patchBoundaryData_[patchName][field] = std::move(bcData);
    }

    /// Set fixed scalar value boundary condition
    void setFixedValue
    (
        const Name& patchName,
        Field field,
        Scalar value
    );

    /// Set fixed scalar gradient boundary condition
    void setFixedGradient
    (
        const Name& patchName,
        Field field,
        Scalar gradient
    );

    /// Set zero gradient boundary condition
    void setZeroGradient
    (
        const Name& patchName,
        Field field
    );

    /// Set no-slip boundary condition
    void setNoSlip
    (
        const Name& patchName,
        Field field
    );

    /// Set wall function boundary condition type
    void setWallFunctionType
    (
        const Name& patchName,
        Field field,
        BCType wallType
    );

// ***************************** Accessor Methods *****************************

    /// Get boundary patch by name
    [[nodiscard]] const BoundaryPatch& patch(const Name& patchName) const;

    /// Get all boundary patches
    [[nodiscard]] const PatchList& patches() const noexcept
    {
        return patches_;
    }

    /// Get number of patches
    [[nodiscard]] Count numPatches() const noexcept
    {
        return patches_.size();
    }

    /// Get boundary condition for a field on a patch
    [[nodiscard]] const BoundaryData& fieldBC
    (
        const Name& patchName,
        Field field
    ) const;

    /// Calculate boundary face value for scalar field
    [[nodiscard]] Scalar boundaryFaceValue
    (
        Field field,
        const ScalarField& phi,
        const Face& boundaryFace
    ) const;

    /// Link boundary faces to their owning patches
    void linkFaces(MutableFaceListRef faces);

    /// Validate boundary condition patch names against mesh patch names
    void validatePatchNames() const;

    /// Print summary of all boundary conditions
    void printSummary() const;

// ****************************** Private Members *****************************

private:

    /// Nested map: patch name → field → boundary data
    PatchBoundaryDataMap patchBoundaryData_;

    /// All boundary patches
    PatchList patches_;

    /// True after linkFaces() — prevents addPatch() after linking
    bool linked_ = false;
};
