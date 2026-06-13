/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file StringTypes.h
 * @brief Intent-revealing owned and borrowed string aliases
 *
 * @details Defines aliases over std::string (owned text) and std::string_view
 * (borrowed text) that signal what the text represents. Borrowed views carry
 * a *Ref suffix that keeps their non-owning nature visible in the name.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <string>
#include <string_view>

// ********************************** Aliases *********************************

/// A name identifier
using Name = std::string;

/// A parser token read from a case or mesh file
using Token = std::string;

/// A path or filename
using FilePath = std::string;

/// Human-readable diagnostic or display text
using Message = std::string;

// ****************************** Borrowed views ******************************

/// A non-owning view of a Name
using NameRef = std::string_view;

/// A non-owning view of a parser Token
using TokenRef = std::string_view;

/// A non-owning view of a FilePath or filename
using FilePathRef = std::string_view;

/// A non-owning view of human-readable diagnostic
using MessageRef = std::string_view;
