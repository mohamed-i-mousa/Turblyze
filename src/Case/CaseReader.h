/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CaseReader.h
 * @brief Parser for the case file
 *
 * @details This header defines the CaseReader class, which parses the input
 * case file for solver configuration.
 *
 * @class CaseReader
 *
 * The CaseReader class provides:
 * - Parsing of the input file
 * - Support for scalar, vector, string, and boolean values
 * - Hierarchical access via dot notation (e.g., "SIMPLE.nIterations")
 * - Automatic type conversion with error handling
 * - Comment and whitespace handling
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <algorithm>
#include <charconv>
#include <cctype>
#include <concepts>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "ErrorHandler.h"
#include "Integer.h"
#include "StringTypes.h"

// *************************** concept CaseInputType **************************

template<typename T>
concept CaseInputType =
    std::same_as<T, Scalar>
 || std::same_as<T, Vector>
 || std::same_as<T, Count>
 || std::same_as<T, Token>
 || std::same_as<T, int>
 || std::same_as<T, bool>;

// ***************************** class CaseReader *****************************

class CaseReader
{
public:

    using NameList = std::vector<Name>;
    using EntryMap = std::map<Name, Token>;
    using SectionMap = std::map<Name, CaseReader>;

// ************************* Special Member Functions *************************

    /// Construct reader from file
    explicit CaseReader(const FilePath& filename);

    /// Default constructor for sections
    CaseReader() noexcept = default;

// ****************************** Public Methods ******************************

    /// Look up a required value
    template<CaseInputType T>
    T lookup(const Name& keyword) const
    {
        const auto it = entries_.find(keyword);
        if (it == entries_.end())
        {
            FatalError("Keyword '" + keyword + "' not found in case file");
        }

        return convertTo<T>(it->second);
    }

    /// Look up an optional value with default
    template<CaseInputType T>
    T lookupOrDefault(const Name& keyword, const T& defaultValue) const
    {
        const auto it = entries_.find(keyword);
        if (it == entries_.end())
        {
            return defaultValue;
        }

        return convertTo<T>(it->second);
    }

    /// Access a section
    const CaseReader& section(const Name& sectionName) const;

    /// Check if section exists
    bool hasSection(const Name& sectionName) const noexcept
    {
        return sections_.contains(sectionName);
    }

    /// Get all section names
    NameList sectionNames() const;

    /// Print case file contents for debugging
    void print(Count indent = 0) const;

// ****************************** Private Methods *****************************

private:

    /// Parse file contents
    void parseFile(const FilePath& filename);

    /// Parse a section from stream
    void parseSection
    (
        std::istream& is,
        CaseReader& sec,
        char terminator = '\0'
    );

    /// Skip comments and whitespace
    void skipCommentsAndWhitespace(std::istream& is);

    /// Read next token from stream
    Token readToken(std::istream& is);

    /// Parse a value as a string representation
    Token parseValue(std::istream& is, const Token& firstToken);

    /// Convert string to specified type
    template<CaseInputType T>
    T convertTo(const Token& value) const
    {
        if constexpr (std::same_as<T, Scalar>)
        {
            Scalar result = S(0.0);
            const auto [ptr, ec] = std::from_chars
            (
                value.data(),
                value.data() + value.size(),
                result
            );

            if (ec != std::errc() || ptr != value.data() + value.size())
            {
                FatalError("Cannot convert '" + value + "' to Scalar");
            }

            return result;
        }
        else if constexpr (std::same_as<T, int>)
        {
            int result = 0;
            const auto [ptr, ec] = std::from_chars
            (
                value.data(),
                value.data() + value.size(),
                result
            );

            if (ec != std::errc() || ptr != value.data() + value.size())
            {
                FatalError("Cannot convert '" + value + "' to int");
            }

            return result;
        }
        else if constexpr (std::same_as<T, Count>)
        {
            Count result = 0;
            const auto [ptr, ec] = std::from_chars
            (
                value.data(),
                value.data() + value.size(),
                result
            );

            if (ec != std::errc() || ptr != value.data() + value.size())
            {
                FatalError("Cannot convert '" + value + "' to Count");
            }

            return result;
        }
        else if constexpr (std::same_as<T, bool>)
        {
            Token lower = value;
            std::ranges::transform(lower, lower.begin(), toLowerSafe);

            if
            (
                lower == "true"
             || lower == "on"
             || lower == "yes"
             || lower == "1"
            )
            {
                return true;
            }

            if
            (
                lower == "false"
             || lower == "off"
             || lower == "no"
             || lower == "0"
            )
            {
                return false;
            }

            FatalError("Cannot convert '" + value + "' to bool");
        }
        else if constexpr (std::same_as<T, Token>)
        {
            return value;
        }
        else if constexpr (std::same_as<T, Vector>)
        {
            // Expecting format: (x y z)
            Token trimmed = value;

            // Remove parentheses
            std::erase(trimmed, '(');
            std::erase(trimmed, ')');

            std::istringstream iss(trimmed);
            Scalar x;
            Scalar y;
            Scalar z;

            if (!(iss >> x >> y >> z))
            {
                FatalError("Cannot convert '" + value + "' to Vector");
            }

            // Reject trailing components: a Vector has exactly three
            Token extra;
            if (iss >> extra)
            {
                FatalError
                (
                    "Cannot convert '" + value
                  + "' to Vector: expected exactly three components"
                );
            }

            return Vector(x, y, z);
        }
    }

    /// Convert a character to lowercase safely
    static char toLowerSafe(char c) noexcept
    {
        return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

// ****************************** Private Members *****************************

private:

    /// Storage for key-value pairs
    EntryMap entries_;

    /// Storage for sections
    SectionMap sections_;

    /// Current file being parsed (for error messages)
    FilePath currentFile_;

    /// Current line number being parsed (for error messages)
    Count currentLine_ = 1;
};
