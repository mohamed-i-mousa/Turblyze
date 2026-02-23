/******************************************************************************
 * @file CaseReader.cpp
 * @brief Implementation of case file parser
 *****************************************************************************/

#include "CaseReader.hpp"
#include <cctype>
#include <iomanip>

// ************************** Constructor/Destructor **************************

CaseReader::CaseReader(const std::string& filename)
{
    currentFile_ = filename;
    parseFile(filename);
}


// ***************************** Accessor Methods *****************************

CaseReader CaseReader::section(const std::string& name) const
{
    auto it = sections_.find(name);
    if (it == sections_.end())
    {
        throw
            std::runtime_error
            (
                "Section '" + name + "' not found"
            );
    }

    return it->second;
}

bool CaseReader::hasSection(const std::string& name) const
{
    return sections_.find(name) != sections_.end();
}

std::vector<std::string> CaseReader::sectionNames() const
{
    std::vector<std::string> names;
    for (const auto& section : sections_)
    {
        names.push_back(section.first);
    }
    return names;
}


// ***************************** Utility Methods ******************************

void CaseReader::print(int indent) const
{
    std::string indentStr(indent * 4, ' ');

    // Print entries
    for (const auto& [key, value] : entries_)
    {
        std::cout
            << indentStr << key << ": " << value << std::endl;
    }

    // Print sections
    for (const auto& [name, section] : sections_)
    {
        std::cout
            << indentStr << name << " {" << std::endl;

        section.print(indent + 1); // recursion

        std::cout
            << indentStr << "}" << std::endl;
    }
}


// ****************************** Private Methods ******************************

void CaseReader::parseFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw
            std::runtime_error
            (
                "Cannot open case file: " + filename
            );
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
    std::string token;

    while ((token = readToken(is)) != "")
    {
        if
        (
            terminator != '\0'
         && token.length() == 1
         && token[0] == terminator
        )
        {
            return;
        }

        // Check for closing brace
        if (token == "}")
        {
            if (terminator == '\0')
            {
                throw
                    std::runtime_error
                    (
                        "Unexpected '}' at line "
                      + std::to_string(currentLine_)
                    );
            }
            return;
        }

        // Read next token to determine if it's a section or key-value pair
        std::string next = readToken(is);

        if (next == "{")
        {
            // It's a section
            CaseReader section;
            parseSection(is, section, '}');
            sec.sections_[token] = section;
        }
        else
        {
            // It's a key-value pair
            std::string value = parseValue(is, next);
            sec.entries_[token] = value;
        }
    }

    // Check if we expected a terminator but hit EOF
    if (terminator != '\0')
    {
        throw
            std::runtime_error
            (
                "Unexpected end of file, expected '"
              + std::string(1, terminator) + "'"
            );
    }
}

void CaseReader::skipCommentsAndWhitespace(std::istream& is) const
{
    char c;
    while (is.get(c))
    {
        if (c == '\n')
        {
            const_cast<int&>(currentLine_)++;
        }

        if (std::isspace(c))
        {
            continue;
        }

        if (c == '/')
        {
            char next = is.peek();
            if (next == '/')
            {
                // Single-line comment
                is.get(); // consume second '/'
                std::string line;
                std::getline(is, line);
                const_cast<int&>(currentLine_)++;
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
                        const_cast<int&>(currentLine_)++;
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

std::string CaseReader::readToken(std::istream& is) const
{
    skipCommentsAndWhitespace(is);

    std::string token;
    char c;

    // Check for special single-character tokens
    if
    (
        is.peek() == '{'
     || is.peek() == '}'
     || is.peek() == ';'
     || is.peek() == '('
    )
    {
        is.get(c);
        return std::string(1, c);
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
                const_cast<int&>(currentLine_)++;
            }
            if (!std::isspace(c))
            {
                is.putback(c);
            }
            break;
        }
        token += c;
    }

    return token;
}

std::string CaseReader::parseValue
(
    std::istream& is,
    const std::string& firstToken
) const
{
    std::string value = firstToken;

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
                const_cast<int&>(currentLine_)++;
            }
            value += c;
            if (c == '(') parenDepth++;
            if (c == ')') parenDepth--;
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
                const_cast<int&>(currentLine_)++;
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
