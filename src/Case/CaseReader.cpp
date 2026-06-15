/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CaseReader.cpp
 * @brief Implementation of case file parser
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "CaseReader.h"

// Standard library headers
#include <fstream>
#include <iostream>

// Project headers
#include "ErrorHandler.h"

// ************************* Special Member Functions *************************

CaseReader::CaseReader(const FilePath& filename)
:
    currentFile_{filename}
{
    parseFile(filename);
}

// ***************************** Accessor Methods *****************************

const CaseReader& CaseReader::section(const Name& sectionName) const
{
    const auto it = sections_.find(sectionName);
    if (it == sections_.end())
    {
        FatalError("Section '" + sectionName + "' not found");
    }

    return it->second;
}


NameList CaseReader::sectionNames() const
{
    NameList sectionList;
    sectionList.reserve(sections_.size());
    for (const auto& section : sections_)
    {
        sectionList.push_back(section.first);
    }
    return sectionList;
}

// ****************************** Public Methods ******************************

void CaseReader::print(Count indent) const
{
    const Message indentStr(indent * 4, ' ');

    // Print entries
    for (const auto& [key, value] : entries_)
    {
        std::cout
            << indentStr << key << ": " << value << '\n';
    }

    // Print sections
    for (const auto& [sectionName, section] : sections_)
    {
        std::cout
            << indentStr << sectionName << " {" << '\n';

        section.print(indent + 1);

        std::cout
            << indentStr << '}' << '\n';
    }
}

// ****************************** Private Methods *****************************

void CaseReader::parseFile(const FilePath& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        FatalError("Cannot open case file: " + filename);
    }

    currentLine_ = 1;
    parseSection(file, *this, '\0');
}


void CaseReader::parseSection
(
    std::istream& is,
    CaseReader& sec,
    char terminator
)
{
    Token t;

    while ((t = readToken(is)) != "")
    {
        if
        (
            terminator != '\0'
         && t.length() == 1
         && t[0] == terminator
        )
        {
            return;
        }

        // Check for closing brace
        if (t == "}")
        {
            if (terminator == '\0')
            {
                FatalError
                (
                    "Unexpected '}' at line "
                  + std::to_string(currentLine_)
                );
            }
            return;
        }

        // Read next token to determine if it's a section or key-value pair
        const Token next = readToken(is);

        if (next == "{")
        {
            // It's a section
            auto& section = sec.sections_[t];
            parseSection(is, section, '}');
        }
        else
        {
            // It's a key-value pair
            Token value = parseValue(is, next);
            sec.entries_[t] = std::move(value);
        }
    }

    // Check if we expected a terminator but hit EOF
    if (terminator != '\0')
    {
        FatalError
        (
            "Unexpected end of file, expected '"
          + Message(1, terminator) + "'"
        );
    }
}


void CaseReader::skipCommentsAndWhitespace(std::istream& is)
{
    char c;
    while (is.get(c))
    {
        if (c == '\n')
        {
            currentLine_++;
        }

        if (std::isspace(c))
        {
            continue;
        }

        if (c == '/')
        {
            const auto next = is.peek();
            if (next == '/')
            {
                // Single-line comment
                is.get(); // consume second '/'
                Token line;
                std::getline(is, line);
                currentLine_++;
                continue;
            }
            else if (next == '*')
            {
                // Multi-line comment
                is.get(); // consume '*'
                char prev = '\0';
                while (is.get(c))
                {
                    if (c == '\n')
                    {
                        currentLine_++;
                    }
                    if (prev == '*' && c == '/')
                    {
                        break;
                    }
                    prev = c;
                }
                continue;
            }
        }

        // Not whitespace or comment, put character back
        is.putback(c);
        break;
    }
}


Token CaseReader::readToken(std::istream& is)
{
    skipCommentsAndWhitespace(is);

    Token t;
    char c;

    // Check for special single-character tokens
    const int next = is.peek();
    if
    (
        next == '{'
     || next == '}'
     || next == ';'
     || next == '('
     || next == ')'
    )
    {
        is.get(c);
        return Token(1, c);
    }

    // Read regular token
    while (is.get(c))
    {
        if
        (
            std::isspace(c)
         || c == '{'
         || c == '}'
         || c == ';'
         || c == '('
         || c == ')'
        )
        {
            if (c == '\n')
            {
                currentLine_++;
            }
            if (!std::isspace(c))
            {
                is.putback(c);
            }
            break;
        }
        t += c;
    }

    return t;
}


Token CaseReader::parseValue
(
    std::istream& is,
    const Token& firstToken
)
{
    Token value = firstToken;

    // Check if it's a vector/list (starts with '(')
    if (firstToken == "(")
    {
        value = "(";
        char c;
        int parenDepth = 1;
        while (is.get(c) && parenDepth > 0)
        {
            if (c == '\n')
            {
                currentLine_++;
            }

            value += c;

            if (c == '(')
            {
                parenDepth++;
            }
            if (c == ')')
            {
                parenDepth--;
            }
        }

        if (parenDepth != 0)
        {
            FatalError
            (
                "Unbalanced parentheses in value '" + value
              + "' in file " + currentFile_
            );
        }

        skipCommentsAndWhitespace(is);
        if (is.peek() == ';')
        {
            is.get(); // consume semicolon
        }
    }
    else
    {
        // Read until semicolon
        skipCommentsAndWhitespace(is);
        char c;
        while (is.get(c))
        {
            if (c == '\n')
            {
                currentLine_++;
            }
            if (c == ';')
            {
                break;
            }
            if (!std::isspace(c))
            {
                value += c;
            }
            else if (!value.empty() && value.back() != ' ')
            {
                value += ' ';
            }
        }
    }

    // Trim trailing whitespace
    while (!value.empty() && std::isspace(value.back()))
    {
        value.pop_back();
    }

    return value;
}
