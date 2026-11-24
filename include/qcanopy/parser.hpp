// qcanopy/parser.hpp — parser for OpenQASM 3.0 and 2.0 circuit files.
//
// Errors are reported as ParseError: 1-based line and column, what() = "line:column: message".
//
// 1. The version header is required, comes first, and selects the dialect:
//
//      OPENQASM 3.0;                       OPENQASM 2.0;
//      include "stdgates.inc";             include "qelib1.inc";
//      qubit[n] name;    (exactly one)     qreg name[n];    (one or more)
//      gphase(a);                          creg name[n];    (ignored)
//      gates: id h x y z s sdg t tdg       gate definitions (point 5)
//             rx ry rz p U cx cz swap      gates: U CX, all of qelib1.inc, the file's own
//
//    U (and 2.0's CX) are built in; other gates need the include first. A circuit has at
//    most kMaxQubits qubits. Register sizes and qubit indices are integer literals.
//    Anything else is an error.
//
// 2. Gates keep source order and operand order (cx: control first). Registers are numbered
//    end to end in declaration order: after "qreg a[2]; qreg b[3];", b[0] is qubit 2.
//
// 3. Gates become IR gates with the same matrix, except:
//    - 3.0's U is e^{i theta/2} times the IR's U, so it adds theta/2 of global phase.
//    - 2.0's rz is read as the IR's RZ, as Qiskit reads it; qelib1.inc's own rz = u1
//      differs from that by a global phase.
//    - 2.0's u1 u2 u3 ry cy ch ccx crz cu1 cu3 are expanded from qelib1.inc's definitions
//      (ccx becomes 15 gates). ch carries a global phase e^{i pi/4}. cu3 has one line added
//      to the specification's text, as in Qiskit's copy, which makes it controlled-U3.
//
// 4. Arguments are expressions over numbers, pi, parentheses, unary '-', and * / + -,
//    evaluated as reals: 1/2 is 0.5 (OpenQASM 3's integer division is not implemented).
//
// 5. A 2.0 definition, "gate name(p, ...) a, b, ... { g(e, ...) a, b; ... }", may call U,
//    CX and gates defined earlier, on its own arguments. Each call is expanded in place.
//    A name can be defined only once.

#pragma once

#include "qcanopy/dag.hpp"
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

        // The same error with the file it came from: what() reads "file:line:column: message".
        ParseError(const std::string& file, const ParseError& e)
            : std::invalid_argument(file + ":" + e.what()), line_(e.line_), column_(e.column_) {}

        uint32_t line() const { return line_; }
        uint32_t column() const { return column_; }

    private:
        uint32_t line_;
        uint32_t column_;
    };

    // Parses circuit source text into a Circuit. Throws ParseError if the text is not
    // accepted; the returned circuit always passes validate().
    Circuit parse_qasm(std::string_view source);

    // As parse_qasm, reading the source from a file. Throws std::runtime_error if the file
    // cannot be read; a ParseError names the file (see ParseError).
    Circuit parse_qasm_file(const std::string& path);

    // The same two entry points, returning the circuit as a DAG for the rewrite passes.
    CircuitDAG parse_qasm_dag(std::string_view source);
    CircuitDAG parse_qasm_file_dag(const std::string& path);

}
