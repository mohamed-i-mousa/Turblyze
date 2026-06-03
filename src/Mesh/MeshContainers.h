/******************************************************************************
 * @file MeshContainers.h
 * @brief Intent-revealing aliases for mesh collections and views
 *
 * @details Defines aliases over std::vector that show mesh elements,
 * and whether it is an owning collection or a borrowed view.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <vector>
#include <span>

// *************************** Forward Declarations ***************************

class Vector;
class Face;
class Cell;
class BoundaryPatch;

// ********************************** Aliases *********************************

/// An ordered list of node coordinates
using NodeList = std::vector<Vector>;

/// An ordered list of faces
using FaceList = std::vector<Face>;

/// An ordered list of cells
using CellList = std::vector<Cell>;

/// An ordered list of boundary patches
using PatchList = std::vector<BoundaryPatch>;

// ****************************** Borrowed views ******************************

/// A non-owning, read-only view of a NodeList
using NodeListRef = std::span<const Vector>;

/// A non-owning, read-only view of a FaceList
using FaceListRef = std::span<const Face>;

/// A non-owning, read-only view of a CellList
using CellListRef = std::span<const Cell>;

/// A non-owning, read-only view of a PatchList
using PatchListRef = std::span<const BoundaryPatch>;

/// A non-owning, mutable view of a FaceList
using MutableFaceListRef = std::span<Face>;

/// A non-owning, mutable view of a CellList
using MutableCellListRef = std::span<Cell>;
