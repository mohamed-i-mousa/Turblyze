/******************************************************************************
 * @file CaseReader.hpp
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

#ifndef CASE_READER_HPP
#define CASE_READER_HPP

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

#include "Scalar.hpp"
#include "Vector.hpp"


class CaseReader
{
public:

    /**
     * @brief Construct reader from file
     * @param filename Path to the input file
     */
    explicit CaseReader(const std::string& filename);

    /// Constructor for sections
    CaseReader() = default;

// Accessor methods

    /**
     * @brief Look up a required value
     * @tparam T Type to convert value
     * @param keyword Key to look up
     * @return Converted value
     * @throws std::runtime_error if keyword not found or conversion fails
     */
    template<typename T>
    T lookup(const std::string& keyword) const;

    /**
     * @brief Look up an optional value with default
     * @tparam T Type to convert value
     * @param keyword Key to look up
     * @param defaultValue Value to return if keyword not found
     * @return Converted value or default
     */
    template<typename T>
    T lookupOrDefault(const std::string& keyword, const T& defaultValue) const;

    /**
     * @brief Access a section
     * @param name Name of section
     * @return Section object
     * @throws std::runtime_error if section not found
     */
    CaseReader section(const std::string& name) const;

    /**
     * @brief Check if section exists
     * @param name Name of section
     * @return true if section exists
     */
    bool hasSection(const std::string& name) const;

    /**
     * @brief Get all section names
     * @return Vector of section names
     */
    std::vector<std::string> sectionNames() const;

    /**
     * @brief Print case file contents (for debugging)
     * @param indent Indentation level
     */
    void print(int indent = 0) const;

private:

// Private methods

    /**
     * @brief Parse file contents
     * @param filename Path to file
     */
    void parseFile(const std::string& filename);

    /**
     * @brief Parse a section from stream
     * @param is Input stream
     * @param sec Section to populate
     * @param terminator Character that ends this section
     */
    void parseSection
    (
        std::istream& is,
        CaseReader& sec,
        char terminator = '\0'
    );

    /**
     * @brief Skip comments and whitespace
     * @param is Input stream
     */
    void skipCommentsAndWhitespace(std::istream& is) const;

    /**
     * @brief Read next token from stream
     * @param is Input stream
     * @return Next token or empty string if EOF
     */
    std::string readToken(std::istream& is) const;

    /**
     * @brief Parse a value (could be scalar, vector, or list)
     * @param is Input stream
     * @param token First token of value
     * @return String representation of value
     */
    std::string parseValue(std::istream& is, const std::string& token) const;

    /**
     * @brief Convert string to specified type
     * @tparam T Target type
     * @param value String value
     * @return Converted value
     */
    template<typename T>
    T convertTo(const std::string& value) const;

// Private members

    /// Storage for key-value pairs
    std::map<std::string, std::string> entries_;

    /// Storage for sections
    std::map<std::string, CaseReader> sections_;

    /// Current file being parsed (for error messages)
    mutable std::string currentFile_;
    mutable int currentLine_ = 1;
};

// Template specializations for common types

template<>
inline Scalar CaseReader::convertTo<Scalar>(const std::string& value) const
{
    try
    {
        return std::stod(value);
    }
    catch (const std::exception& e)
    {
        throw   
            std::runtime_error
            (
                "Cannot convert '" + value + "' to Scalar"
            );
    }
}

template<>
inline int CaseReader::convertTo<int>(const std::string& value) const
{
    try
    {
        return std::stoi(value);
    }
    catch (const std::exception& e)
    {
        throw   
            std::runtime_error
            (
                "Cannot convert '" + value + "' to int"
            );
    }
}

template<>
inline bool CaseReader::convertTo<bool>(const std::string& value) const
{
    std::string lower = value;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    if (lower == "true" || lower == "on" || lower == "yes" || lower == "1")
    {
        return true;
    }
    else if
    (
        lower == "false" || lower == "off" || lower == "no" || lower == "0"
    )
    {
        return false;
    }

    throw   
        std::runtime_error
        (
            "Cannot convert '" + value + "' to bool"
        );
}

template<>
inline std::string CaseReader::convertTo<std::string>
(
    const std::string& value
) const
{
    return value;
}

template<>
inline Vector CaseReader::convertTo<Vector>(const std::string& value) const
{
    // Expecting format: (x y z)
    std::string trimmed = value;

    // Remove parentheses and whitespace
    trimmed.erase
    (
        std::remove(trimmed.begin(), trimmed.end(), '('), trimmed.end()
    );

    trimmed.erase
    (
        std::remove(trimmed.begin(), trimmed.end(), ')'), trimmed.end()
    );

    std::istringstream iss(trimmed);
    Scalar x, y, z;

    if (!(iss >> x >> y >> z))
    {
        throw   
            std::runtime_error
            (
                "Cannot convert '" + value + "' to Vector"
            );
    }

    return Vector(x, y, z);
}

// Template method implementations

template<typename T>
T CaseReader::lookup(const std::string& keyword) const
{
    auto it = entries_.find(keyword);
    if (it == entries_.end())
    {
        throw   
            std::runtime_error
            (
                "Keyword '" + keyword + "' not found in case file"
            );
    }

    return convertTo<T>(it->second);
}

template<typename T>
T CaseReader::lookupOrDefault
(
    const std::string& keyword, const T& defaultValue
) const
{
    auto it = entries_.find(keyword);
    if (it == entries_.end())
    {
        return defaultValue;
    }

    try
    {
        return convertTo<T>(it->second);
    }
    catch (const std::exception&)
    {
        return defaultValue;
    }
}

#endif // CASE_READER_HPP
