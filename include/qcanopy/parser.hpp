// qcanopy/parser.hpp — circuit source parser.
//
// Errors in the source text are reported as ParseError, which carries the 1-based line and
// column where the problem starts. what() reads "line:column: message".

#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

namespace qcanopy {

    class ParseError : public std::invalid_argument {
    public:
        ParseError(uint32_t line, uint32_t column, const std::string& message)
            : std::invalid_argument(std::to_string(line) + ":" + std::to_string(column) + ": " + message),
              line_(line), column_(column) {}

        uint32_t line() const { return line_; }
        uint32_t column() const { return column_; }

    private:
        uint32_t line_;
        uint32_t column_;
    };

}
