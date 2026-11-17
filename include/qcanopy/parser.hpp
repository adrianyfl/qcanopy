// qcanopy/parser.hpp — parser for OpenQASM 3.0 and 2.0 circuit files.
//
// Errors in the source text are reported as ParseError, which carries the 1-based line and
// column where the problem starts. what() reads "line:column: message".
//
// 1. Accepted source is a subset of OpenQASM 3.0 or of OpenQASM 2.0. The version header is
//    required, comes first, and selects which:
//
//      OPENQASM 3.0;                       OPENQASM 2.0;
//      include "stdgates.inc";             include "qelib1.inc";
//      qubit[n] name;    (exactly one)     qreg name[n];    (one or more)
//      gphase(a);                          creg name[n];    (checked, then ignored)
//      gates: id h x y z s sdg t tdg       gates: id x y z h s sdg t tdg
//             rx ry rz p U cx cz swap             rx rz cx cz
//
//    The include must precede any gate other than 3.0's U. A gate statement is
//      g q[i];   g(a) q[i];   g q[i],q[j];
//    in both versions, and gphase(a) adds a to Circuit::global_phase. The registers hold
//    1 to kMaxQubits qubits in total. Arguments are constant expressions (point 4);
//    register sizes and qubit indices are integer literals. Comments, "// ..." and
//    "/* ... */", may appear wherever whitespace may.
//
// 2. Gates are added in source order, and operand order is kept: for cx the first operand
//    is the control (ir.hpp, point 3). Registers are laid end to end in declaration order:
//    after "qreg a[2]; qreg b[3];" the qubit b[0] is qubit 2 of the Circuit.
//
// 3. Every gate maps to the IR gate with the same matrix, with two exceptions.
//    3.0's built-in U(theta, phi, lambda) is e^{i theta/2} times the IR's phase-free U, so
//    it becomes a U gate plus theta/2 of global phase.
//    2.0's qelib1.inc defines rz(phi) as u1(phi), which is e^{i phi/2} times the IR's RZ.
//    It is read as plain RZ with no phase added, as Qiskit reads it, so a 2.0 circuit
//    with rz gates is exact only up to a global phase.
//
// 4. A constant expression is built from numeric literals, the constant pi, parentheses,
//    unary '-', and the operators * / and + - with their usual precedence. Everything is
//    evaluated as a real, so 1/2 is 0.5. (OpenQASM 3 divides two integers as integers,
//    giving 0; that rule is not implemented.)

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
