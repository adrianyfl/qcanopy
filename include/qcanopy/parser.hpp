// qcanopy/parser.hpp — circuit source parser.
//
// Errors in the source text are reported as ParseError, which carries the 1-based line and
// column where the problem starts. what() reads "line:column: message".
//
// 1. Accepted source, in this order:
//      version header, 3.0                 required, first
//      include "stdgates.inc";             required before any gate other than U
//      qubit[n] name;                      exactly one register, 1 <= n <= kMaxQubits
//    followed by any number of
//      gphase(a);                          adds a to Circuit::global_phase
//      g q[i];   g(a) q[i];   g q[i],q[j];
//    where g is one of id h x y z s sdg t tdg rx ry rz p U cx cz swap.
//    Arguments are numeric literals with an optional leading '-'.
//
// 2. Gates are added in source order, and operand order is kept: for cx the first operand
//    is the control (ir.hpp, point 3).
//
// 3. Every gate maps to the IR gate with the same matrix except U: the source's
//    U(theta, phi, lambda) is e^{i theta/2} times the IR's phase-free U, so it becomes a
//    U gate plus theta/2 of global phase.

#pragma once

#include "qcanopy/ir.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

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

    // Parses circuit source text into a Circuit. Throws ParseError if the text is not
    // accepted; the returned circuit always passes validate().
    Circuit parse_qasm(std::string_view source);

}
