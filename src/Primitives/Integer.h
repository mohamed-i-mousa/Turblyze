/******************************************************************************
 * @file Integer.h
 * @brief Intent-revealing integer aliases
 *
 * @details Defines two aliases over std::size_t that signal how an integer is
 * used: Index addresses an element in storage, Count is a size or quantity.
 * The distinction is purely for the reader.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <cstddef>
#include <vector>
#include <span>

// ********************************** Aliases *********************************

/// An entity or array index
using Index = std::size_t;

/// A size or quantity
using Count = std::size_t;

/// An ordered, list of indices
using IndexList = std::vector<Index>;

/// An ordered, list of counts
using CountList = std::vector<Count>;

// ****************************** Borrowed views ******************************

/// A non-owning, read-only view of an IndexList
using IndexListRef = std::span<const Index>;
