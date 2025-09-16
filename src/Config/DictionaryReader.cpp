/******************************************************************************
 * @file DictionaryReader.cpp
 * @brief Implementation of OpenFOAM-style dictionary parser
 *
 * @author Mohamed Mousa
 * @date 2025
 *****************************************************************************/

#include "DictionaryReader.hpp"
#include <cctype>
#include <iomanip>

DictionaryReader::DictionaryReader(const std::string& filename)
{
    currentFile_ = filename;
    parseFile(filename);
}

void DictionaryReader::parseFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open dictionary file: " + filename);
    }

    currentLine_ = 1;
    parseDict(file, *this, '\0');
}

void DictionaryReader::parseDict(std::istream& is, DictionaryReader& dict, char terminator)
{
    std::string token;

    while ((token = readToken(is)) != "") {
        // Check for dictionary terminator
        if (terminator != '\0' && token.length() == 1 && token[0] == terminator) {
            return;
        }

        // Check for closing brace (shouldn't happen at top level)
        if (token == "}") {
            if (terminator == '\0') {
                throw std::runtime_error("Unexpected '}' at line " +
                                       std::to_string(currentLine_));
            }
            return;
        }

        // Read next token to determine if it's a sub-dictionary or key-value pair
        std::string next = readToken(is);

        if (next == "{") {
            // It's a sub-dictionary
            DictionaryReader subDict;
            parseDict(is, subDict, '}');
            dict.subDicts_[token] = subDict;
        } else {
            // It's a key-value pair
            std::string value = parseValue(is, next);
            dict.entries_[token] = value;
        }
    }

    // Check if we expected a terminator but hit EOF
    if (terminator != '\0') {
        throw std::runtime_error("Unexpected end of file, expected '" +
                               std::string(1, terminator) + "'");
    }
}

void DictionaryReader::skipCommentsAndWhitespace(std::istream& is) const
{
    char c;
    while (is.get(c)) {
        if (c == '\n') {
            const_cast<int&>(currentLine_)++;
        }

        if (std::isspace(c)) {
            continue;
        }

        if (c == '/') {
            char next = is.peek();
            if (next == '/') {
                // Single-line comment
                is.get(); // consume second '/'
                std::string line;
                std::getline(is, line);
                const_cast<int&>(currentLine_)++;
                continue;
            } else if (next == '*') {
                // Multi-line comment
                is.get(); // consume '*'
                char prev = '\0';
                while (is.get(c)) {
                    if (c == '\n') {
                        const_cast<int&>(currentLine_)++;
                    }
                    if (prev == '*' && c == '/') {
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

std::string DictionaryReader::readToken(std::istream& is) const
{
    skipCommentsAndWhitespace(is);

    std::string token;
    char c;

    // Check for special single-character tokens
    if (is.peek() == '{' || is.peek() == '}' || is.peek() == ';' || is.peek() == '(') {
        is.get(c);
        return std::string(1, c);
    }

    // Read regular token
    while (is.get(c)) {
        if (std::isspace(c) || c == '{' || c == '}' || c == ';' || c == '(' || c == ')') {
            if (c == '\n') {
                const_cast<int&>(currentLine_)++;
            }
            if (!std::isspace(c)) {
                is.putback(c);
            }
            break;
        }
        token += c;
    }

    return token;
}

std::string DictionaryReader::parseValue(std::istream& is, const std::string& firstToken) const
{
    std::string value = firstToken;

    // Check if it's a vector/list (starts with '(')
    if (firstToken == "(") {
        value = "(";
        char c;
        int parenDepth = 1;
        while (is.get(c) && parenDepth > 0) {
            if (c == '\n') {
                const_cast<int&>(currentLine_)++;
            }
            value += c;
            if (c == '(') parenDepth++;
            if (c == ')') parenDepth--;
        }
    } else {
        // Read until semicolon
        skipCommentsAndWhitespace(is);
        char c;
        while (is.get(c)) {
            if (c == '\n') {
                const_cast<int&>(currentLine_)++;
            }
            if (c == ';') {
                break;
            }
            if (!std::isspace(c)) {
                value += c;
            } else if (!value.empty() && value.back() != ' ') {
                value += ' ';
            }
        }
    }

    // Trim trailing whitespace
    while (!value.empty() && std::isspace(value.back())) {
        value.pop_back();
    }

    return value;
}

DictionaryReader DictionaryReader::subDict(const std::string& name) const
{
    auto it = subDicts_.find(name);
    if (it == subDicts_.end()) {
        throw std::runtime_error("Sub-dictionary '" + name + "' not found");
    }
    return it->second;
}

bool DictionaryReader::found(const std::string& keyword) const
{
    return entries_.find(keyword) != entries_.end();
}

bool DictionaryReader::foundSubDict(const std::string& name) const
{
    return subDicts_.find(name) != subDicts_.end();
}

std::vector<std::string> DictionaryReader::keywords() const
{
    std::vector<std::string> keys;
    for (const auto& entry : entries_) {
        keys.push_back(entry.first);
    }
    return keys;
}

std::vector<std::string> DictionaryReader::subDictNames() const
{
    std::vector<std::string> names;
    for (const auto& subDict : subDicts_) {
        names.push_back(subDict.first);
    }
    return names;
}

void DictionaryReader::print(int indent) const
{
    std::string indentStr(indent * 4, ' ');

    // Print entries
    for (const auto& [key, value] : entries_) {
        std::cout << indentStr << key << ": " << value << std::endl;
    }

    // Print sub-dictionaries
    for (const auto& [name, subDict] : subDicts_) {
        std::cout << indentStr << name << " {" << std::endl;
        subDict.print(indent + 1);
        std::cout << indentStr << "}" << std::endl;
    }
}