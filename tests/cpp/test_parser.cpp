// Tests for parse_qasm (qcanopy/parser.hpp).
#include "qcanopy/generator.hpp"
#include "qcanopy/parser.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using namespace qcanopy;

namespace {
    int failures = 0;

#define CHECK(cond) \
    do { if (!(cond)) { std::printf("%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)

    using K = GateKind;

    // Three lines, so a body appended to it starts on line 4.
    const std::string kHeader = "OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[3] q;\n";

    Circuit parse_body(const std::string& body) { return parse_qasm(kHeader + body); }

    // Same kind, qubits (in order), matrix and parameters.
    bool is_gate(const Gate& got, const Gate& want) {
        return same_gate(got, want, 0.0) && got.q == want.q && got.params == want.params;
    }

    // True if parsing src throws a ParseError positioned at (line, column).
    bool fails_at(const std::string& src, uint32_t line, uint32_t column) {
        try {
            (void)parse_qasm(src);
        } catch (const ParseError& e) {
            return e.line() == line && e.column() == column;
        }
        return false;
    }

    // As fails_at, for a body placed after kHeader.
    bool body_fails_at(const std::string& body, uint32_t line, uint32_t column) {
        return fails_at(kHeader + body, line, column);
    }

    std::string error_of(const std::string& src) {
        try {
            (void)parse_qasm(src);
        } catch (const std::invalid_argument& e) {      // ParseError is an invalid_argument
            return e.what();
        }
        return "";
    }

    // A program in the layout qcanopy_generate --format qasm writes.
    void test_program() {
        const Circuit c = parse_qasm(
            "OPENQASM 3.0;\n"
            "include \"stdgates.inc\";\n"
            "qubit[4] q;\n"
            "gphase(4.3196898986859651);\n"
            "rz(-0.78539816339744828) q[0];\n"
            "rx(-1.5707963267948966) q[0];\n"
            "s q[0];\n"
            "cx q[1],q[0];\n"
            "swap q[2],q[3];\n");
        validate(c);
        CHECK(c.n == 4);
        CHECK(c.global_phase == 4.3196898986859651);
        CHECK(c.size() == 5);
        if (c.size() != 5) return;
        CHECK(is_gate(c.gates[0], make_gate(K::RZ, 0, {-0.78539816339744828, 0, 0})));
        CHECK(is_gate(c.gates[1], make_gate(K::RX, 0, {-1.5707963267948966, 0, 0})));
        CHECK(is_gate(c.gates[2], make_gate(K::S, 0)));
        CHECK(is_gate(c.gates[3], make_gate(K::CX, 1, 0)));
        CHECK(is_gate(c.gates[4], make_gate(K::SWAP, 2, 3)));
    }

    void test_every_gate() {
        struct Case {
            const char* source;
            Gate want;
        };
        const std::vector<Case> cases = {
            {"id q[1];",                make_gate(K::I, 1)},
            {"h q[1];",                 make_gate(K::H, 1)},
            {"x q[1];",                 make_gate(K::X, 1)},
            {"y q[1];",                 make_gate(K::Y, 1)},
            {"z q[1];",                 make_gate(K::Z, 1)},
            {"s q[1];",                 make_gate(K::S, 1)},
            {"sdg q[1];",               make_gate(K::Sdg, 1)},
            {"t q[1];",                 make_gate(K::T, 1)},
            {"tdg q[1];",               make_gate(K::Tdg, 1)},
            {"rx(0.25) q[1];",          make_gate(K::RX, 1, {0.25, 0, 0})},
            {"ry(0.25) q[1];",          make_gate(K::RY, 1, {0.25, 0, 0})},
            {"rz(0.25) q[1];",          make_gate(K::RZ, 1, {0.25, 0, 0})},
            {"p(0.25) q[1];",           make_gate(K::P, 1, {0.25, 0, 0})},
            {"U(0.25,0.5,0.75) q[1];",  make_gate(K::U, 1, {0.25, 0.5, 0.75})},
            {"cx q[2],q[0];",           make_gate(K::CX, 2, 0)},
            {"cz q[2],q[0];",           make_gate(K::CZ, 2, 0)},
            {"swap q[2],q[0];",         make_gate(K::SWAP, 2, 0)},
        };
        for (const Case& k : cases) {
            const Circuit c = parse_body(k.source);
            validate(c);
            const bool ok = c.n == 3 && c.size() == 1 && is_gate(c.gates[0], k.want);
            if (!ok) std::printf("gate case failed: %s\n", k.source);
            CHECK(ok);
        }
    }

    // parser.hpp, point 3: U carries theta/2 of global phase; no other gate carries any.
    void test_global_phase() {
        CHECK(parse_body("").global_phase == 0.0);
        CHECK(parse_body("h q[0]; rz(0.5) q[1]; rx(0.5) q[2]; p(0.5) q[0]; cx q[0],q[1];").global_phase == 0.0);

        CHECK(parse_body("U(0.5,0.2,0.3) q[0];").global_phase == 0.25);
        CHECK(parse_body("U(0.5,0.2,0.3) q[0]; U(-1,0,0) q[1];").global_phase == -0.25);
        // The writer emits gphase(-theta/2) for an IR U gate; reading it back cancels.
        CHECK(parse_body("gphase(-0.25); U(0.5,0.2,0.3) q[0];").global_phase == 0.0);

        CHECK(parse_body("gphase(1.5); gphase(-0.5); h q[0]; gphase(2);").global_phase == 3.0);
        // gphase needs no qubits, so it may precede the register.
        CHECK(parse_qasm("OPENQASM 3.0;\ngphase(0.5);\nqubit[1] q;\n").global_phase == 0.5);
    }

    // The value an expression has as a gate argument.
    double angle(const std::string& text) {
        const Circuit c = parse_body("rz(" + text + ") q[0];");
        return c.size() == 1 ? c.gates[0].params[0] : -999.0;
    }

    constexpr double kPi = 3.14159265358979323846;

    // parser.hpp, point 4.
    void test_expressions() {
        CHECK(angle("pi") == kPi);
        CHECK(angle("pi/2") == kPi / 2);
        CHECK(angle("-pi/4") == -kPi / 4);
        CHECK(angle("3*pi") == 3 * kPi);
        CHECK(angle("-3*pi/2") == -3 * kPi / 2);
        CHECK(angle("0.5*pi") == 0.5 * kPi);
        CHECK(angle("6.103515625e-05*pi") == 6.103515625e-05 * kPi);
        CHECK(angle("pi - pi/8 + 1") == kPi - kPi / 8 + 1);

        // Precedence and associativity.
        CHECK(angle("1+2*3") == 7.0);
        CHECK(angle("(1+2)*3") == 9.0);
        CHECK(angle("10-4-3") == 3.0);
        CHECK(angle("2/4/2") == 0.25);
        CHECK(angle("2*-3") == -6.0);
        CHECK(angle("--1") == 1.0);
        CHECK(angle("-(-(pi))") == kPi);
        CHECK(angle("((((1))))") == 1.0);

        // Integer literals are reals like any other number.
        CHECK(angle("1/2") == 0.5);
        CHECK(angle("-7/2") == -3.5);
        CHECK(angle("3/2*pi") == 1.5 * kPi);
        CHECK(angle("3*pi/2") == 3 * kPi / 2);

        // The same expressions work for gphase and for every parameter of U.
        CHECK(parse_body("gphase(pi/4); gphase(-pi/4 + 1);").global_phase == kPi / 4 + (-kPi / 4 + 1));
        const Circuit u = parse_body("U(pi/2, 0, 2*pi) q[1];");
        CHECK(u.size() == 1 && is_gate(u.gates[0], make_gate(K::U, 1, {kPi / 2, 0, 2 * kPi})));
        CHECK(u.global_phase == kPi / 4);
    }

    void test_expression_errors() {
        // A result that is not a finite number is reported where the expression starts.
        CHECK(body_fails_at("rz(1/0) q[0];", 4, 4));
        CHECK(body_fails_at("rz(pi/(1-1)) q[0];", 4, 4));
        CHECK(body_fails_at("rz(0/0) q[0];", 4, 4));
        CHECK(body_fails_at("rz(1e300*1e300) q[0];", 4, 4));

        CHECK(body_fails_at("rz(foo) q[0];", 4, 4));
        CHECK(body_fails_at("rz(q) q[0];", 4, 4));                  // a register is not a number
        CHECK(body_fails_at("rz(tau) q[0];", 4, 4));                // pi is the only constant
        CHECK(body_fails_at("rz(sin(1)) q[0];", 4, 4));             // no functions
        CHECK(body_fails_at("rz(2**3) q[0];", 4, 5));               // no power operator
        CHECK(body_fails_at("rz(+1) q[0];", 4, 4));                 // no unary plus
        CHECK(body_fails_at("rz(*2) q[0];", 4, 4));
        CHECK(body_fails_at("rz(1+) q[0];", 4, 6));
        CHECK(body_fails_at("rz(1 2) q[0];", 4, 6));
        CHECK(body_fails_at("rz((1) q[0];", 4, 8));                 // unbalanced '('
        CHECK(body_fails_at("rz(1)) q[0];", 4, 6));
        CHECK(body_fails_at("rz(2^3) q[0];", 4, 5));

        // Register sizes and qubit indices are literals, not expressions.
        CHECK(fails_at("OPENQASM 3.0;\nqubit[1+1] q;", 2, 8));
        CHECK(body_fails_at("h q[1-1];", 4, 6));

        CHECK(error_of(kHeader + "rz(foo) q[0];") == "4:4: unknown constant 'foo'");
        CHECK(error_of(kHeader + "rz(1/0) q[0];") ==
              "4:4: expression does not evaluate to a finite number");
    }

    void test_numbers() {
        CHECK(angle("1") == 1.0);
        CHECK(angle("-1") == -1.0);
        CHECK(angle("0") == 0.0);
        CHECK(angle("5.") == 5.0);
        CHECK(angle(".5") == 0.5);
        CHECK(angle("-.5") == -0.5);
        CHECK(angle("1e-05") == 1e-05);
        CHECK(angle("-2.5E+2") == -250.0);
        CHECK(angle("0.78539816339744828") == 0.78539816339744828);
        CHECK(angle("  - 3.0 ") == -3.0);
    }

    void test_layout() {
        // Statements need not be one per line; whitespace between tokens is free.
        const Circuit a = parse_qasm("OPENQASM 3.0;include\"stdgates.inc\";qubit[2]q;h q[0];cx q[0],q[1];");
        CHECK(a.n == 2 && a.size() == 2);

        const Circuit b = parse_qasm("\r\n  OPENQASM   3.0 ;\r\n\tinclude \"stdgates.inc\" ;\r\n"
                                     "qubit [ 2 ] q ;\r\ncx\tq [ 0 ] , q [ 1 ] ;\r\n");
        CHECK(b.n == 2 && b.size() == 1 && is_gate(b.gates[0], make_gate(K::CX, 0, 1)));

        // "3" is the same version as "3.0"; the register may have any name; empty
        // parentheses are allowed on a gate without parameters.
        const Circuit c = parse_qasm("OPENQASM 3;\ninclude \"stdgates.inc\";\nqubit[1] wires;\nh() wires[0];\n");
        CHECK(c.n == 1 && c.size() == 1 && is_gate(c.gates[0], make_gate(K::H, 0)));

        // U is part of the language, so it needs no include.
        const Circuit d = parse_qasm("OPENQASM 3.0;\nqubit[1] q;\nU(0,0,0) q[0];\n");
        CHECK(d.size() == 1 && d.gates[0].kind == K::U);

        // A circuit may have no gates, and may use every qubit the IR allows.
        const Circuit e = parse_qasm("OPENQASM 3.0;\nqubit[32] q;\nU(0,0,0) q[31];\n");
        validate(e);
        CHECK(e.n == 32 && e.gates[0].q[0] == 31);
        CHECK(parse_body("").size() == 0);
    }

    void test_comments() {
        const Circuit c = parse_qasm(
            "// generated file\n"
            "OPENQASM 3.0; // version\n"
            "include \"stdgates.inc\";\n"
            "/* two qubits,\n"
            "   one register */\n"
            "qubit[2] q;\n"
            "h q[0]; // superposition\n"
            "cx /* control */ q[0], /* target */ q[1];\n"
            "// rz(0.25) q[0];\n"
            "rz(/* angle */ 0.5) q[1]; /* trailing */");
        CHECK(c.n == 2 && c.size() == 3);
        if (c.size() == 3) {
            CHECK(is_gate(c.gates[0], make_gate(K::H, 0)));
            CHECK(is_gate(c.gates[1], make_gate(K::CX, 0, 1)));
            CHECK(is_gate(c.gates[2], make_gate(K::RZ, 1, {0.5, 0, 0})));
        }

        // Commented-out statements are not parsed, so they cannot fail.
        CHECK(parse_body("// h q[99];\n/* foo q[0]; */").size() == 0);

        // Errors after a comment keep their real position.
        CHECK(body_fails_at("/* a\n b */ h q[9];", 5, 11));
        CHECK(body_fails_at("h q[0]; /* never closed", 4, 9));
    }

    void test_dag() {
        const CircuitDAG dag = parse_qasm_dag(kHeader + "gphase(0.5); h q[0]; cx q[0],q[1]; t q[2];");
        CHECK(dag.n == 3);
        CHECK(dag.global_phase == 0.5);
        CHECK(dag.nodes.size() == 3);
        if (dag.nodes.size() != 3) return;

        // Nodes are in source order: 0 = h, 1 = cx, 2 = t.
        CHECK(is_gate(dag.nodes[0].gate, make_gate(K::H, 0)));
        CHECK(is_gate(dag.nodes[1].gate, make_gate(K::CX, 0, 1)));
        CHECK(is_gate(dag.nodes[2].gate, make_gate(K::T, 2)));
        // Wire 0 runs h -> cx, wire 1 holds only cx, wire 2 only t.
        CHECK(dag.head[0] == 0 && dag.next_on(0, 0) == 1 && dag.tail[0] == 1);
        CHECK(dag.head[1] == 1 && dag.tail[1] == 1);
        CHECK(dag.head[2] == 2 && dag.tail[2] == 2);

        try {
            (void)parse_qasm_dag(kHeader + "h q[7];");
            CHECK(false);
        } catch (const ParseError& e) {
            CHECK(e.line() == 4 && e.column() == 5);
        }
    }

    void test_header_errors() {
        CHECK(fails_at("", 1, 1));
        CHECK(fails_at("qubit[2] q;", 1, 1));                               // no version
        CHECK(fails_at("OPENQASM 1.0;", 1, 10));
        CHECK(fails_at("OPENQASM 3.1;", 1, 10));
        CHECK(fails_at("OPENQASM;", 1, 9));
        CHECK(fails_at("OPENQASM 3.0\nqubit[2] q;", 2, 1));                 // missing ';'
        CHECK(fails_at("OPENQASM 3.0;\nOPENQASM 3.0;", 2, 1));

        CHECK(fails_at("OPENQASM 3.0;\ninclude \"qelib1.inc\";", 2, 9));
        CHECK(fails_at("OPENQASM 3.0;\ninclude stdgates.inc;", 2, 9));
        CHECK(fails_at("OPENQASM 3.0;\nqubit[1] q;\nh q[0];", 3, 1));       // h before the include

        CHECK(fails_at("OPENQASM 3.0;\ninclude \"stdgates.inc\";\n", 3, 1));    // no register
        CHECK(fails_at("OPENQASM 3.0;\ninclude \"stdgates.inc\";\nh q[0];\nqubit[1] q;", 3, 3));
        CHECK(fails_at("OPENQASM 3.0;\nqubit[0] q;", 2, 7));
        CHECK(fails_at("OPENQASM 3.0;\nqubit[33] q;", 2, 7));
        CHECK(fails_at("OPENQASM 3.0;\nqubit[2.0] q;", 2, 7));
        CHECK(fails_at("OPENQASM 3.0;\nqubit q;", 2, 7));                   // unsized register
        CHECK(fails_at("OPENQASM 3.0;\nqubit[2];", 2, 9));
        CHECK(body_fails_at("qubit[2] r;", 4, 1));                          // second register
    }

    void test_gate_errors() {
        CHECK(body_fails_at("foo q[0];", 4, 1));
        CHECK(body_fails_at("ccx q[0],q[1],q[2];", 4, 1));
        CHECK(body_fails_at(";", 4, 1));
        CHECK(body_fails_at("3 q[0];", 4, 1));

        // Operands.
        CHECK(body_fails_at("h q[3];", 4, 5));                              // register has 3 qubits
        CHECK(body_fails_at("h q[99999999999999999999];", 4, 5));
        CHECK(body_fails_at("h q[-1];", 4, 5));
        CHECK(body_fails_at("h r[0];", 4, 3));
        CHECK(body_fails_at("h q;", 4, 4));                                 // whole register
        CHECK(body_fails_at("h q[0]", 4, 7));                               // missing ';'
        CHECK(body_fails_at("h;", 4, 2));
        CHECK(body_fails_at("cx q[0] q[1];", 4, 9));                        // missing ','
        CHECK(body_fails_at("cx q[1],q[1];", 4, 9));
        CHECK(body_fails_at("cx q[0];", 4, 1));
        CHECK(body_fails_at("h q[0],q[1];", 4, 1));
        CHECK(body_fails_at("cx q[0],q[1],q[2];", 4, 1));

        // Parameters.
        CHECK(body_fails_at("rz q[0];", 4, 1));
        CHECK(body_fails_at("rz() q[0];", 4, 1));
        CHECK(body_fails_at("rz(1,2) q[0];", 4, 1));
        CHECK(body_fails_at("h(1) q[0];", 4, 1));
        CHECK(body_fails_at("U(1,2) q[0];", 4, 1));
        CHECK(body_fails_at("U(1,2,3,4) q[0];", 4, 1));
        CHECK(body_fails_at("rz(1 q[0];", 4, 6));
        CHECK(body_fails_at("rz(1,) q[0];", 4, 6));
        CHECK(body_fails_at("rz(1e999) q[0];", 4, 4));

        CHECK(body_fails_at("gphase 1;", 4, 8));
        CHECK(body_fails_at("gphase(1) q[0];", 4, 11));
        CHECK(body_fails_at("gphase();", 4, 8));

        // Lexer errors come through unchanged.
        CHECK(body_fails_at("h q[0]; @", 4, 9));
        CHECK(body_fails_at("rz(1e) q[0];", 4, 4));
    }

    // ---- OpenQASM 2.0 ----

    // Three lines, so a body appended to it starts on line 4.
    const std::string kHeader2 = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\n";

    // A program in the layout of the benchmark circuits: a leading comment, two quantum
    // registers, an unused classical register, and angles written with pi.
    void test_qasm2_program() {
        const Circuit c = parse_qasm(
            "// Used Gate Set: ['rz', 'rx', 'x', 'cx']\n"
            "OPENQASM 2.0;\n"
            "include \"qelib1.inc\";\n"
            "qreg q[2];\n"
            "qreg flag[1];\n"
            "creg c[2];\n"
            "s q[0];\n"
            "rx(pi/2) q[0];\n"
            "rz(-3*pi/2) q[1];\n"
            "x flag[0];\n"
            "cx q[1],flag[0];\n"
            "cz q[0], q[1];\n");
        validate(c);
        CHECK(c.n == 3);
        CHECK(c.global_phase == 0.0);               // parser.hpp, point 3: rz adds no phase
        CHECK(c.size() == 6);
        if (c.size() != 6) return;
        CHECK(is_gate(c.gates[0], make_gate(K::S, 0)));
        CHECK(is_gate(c.gates[1], make_gate(K::RX, 0, {kPi / 2, 0, 0})));
        CHECK(is_gate(c.gates[2], make_gate(K::RZ, 1, {-3 * kPi / 2, 0, 0})));
        CHECK(is_gate(c.gates[3], make_gate(K::X, 2)));         // flag[0] follows q[0..1]
        CHECK(is_gate(c.gates[4], make_gate(K::CX, 1, 2)));
        CHECK(is_gate(c.gates[5], make_gate(K::CZ, 0, 1)));
    }

    void test_qasm2_gates() {
        struct Case {
            const char* source;
            Gate want;
        };
        const std::vector<Case> cases = {
            {"id q[1];",        make_gate(K::I, 1)},
            {"x q[1];",         make_gate(K::X, 1)},
            {"y q[1];",         make_gate(K::Y, 1)},
            {"z q[1];",         make_gate(K::Z, 1)},
            {"h q[1];",         make_gate(K::H, 1)},
            {"s q[1];",         make_gate(K::S, 1)},
            {"sdg q[1];",       make_gate(K::Sdg, 1)},
            {"t q[1];",         make_gate(K::T, 1)},
            {"tdg q[1];",       make_gate(K::Tdg, 1)},
            {"rx(0.25) q[1];",  make_gate(K::RX, 1, {0.25, 0, 0})},
            {"rz(0.25) q[1];",  make_gate(K::RZ, 1, {0.25, 0, 0})},
            {"cx q[2],q[0];",   make_gate(K::CX, 2, 0)},
            {"cz q[2],q[0];",   make_gate(K::CZ, 2, 0)},
        };
        for (const Case& k : cases) {
            const Circuit c = parse_qasm(kHeader2 + k.source);
            validate(c);
            const bool ok = c.n == 3 && c.size() == 1 && is_gate(c.gates[0], k.want) && c.global_phase == 0.0;
            if (!ok) std::printf("2.0 gate case failed: %s\n", k.source);
            CHECK(ok);
        }
    }

    // parser.hpp, point 2: registers are numbered end to end in declaration order.
    void test_qasm2_registers() {
        const std::string head = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n";

        const Circuit c = parse_qasm(head + "qreg a[2];\nqreg b[3];\nqreg c[1];\n"
                                            "x a[1]; x b[0]; x b[2]; x c[0]; cx c[0],a[0];");
        CHECK(c.n == 6 && c.size() == 5);
        if (c.size() == 5) {
            CHECK(c.gates[0].q[0] == 1);
            CHECK(c.gates[1].q[0] == 2);
            CHECK(c.gates[2].q[0] == 4);
            CHECK(c.gates[3].q[0] == 5);
            CHECK(is_gate(c.gates[4], make_gate(K::CX, 5, 0)));
        }

        // A register may be declared after gates have used earlier ones.
        const Circuit late = parse_qasm(head + "qreg a[1];\nx a[0];\nqreg b[1];\ncx a[0],b[0];");
        CHECK(late.n == 2 && late.size() == 2);

        // Classical registers do not count as qubits, wherever they appear.
        const Circuit classical = parse_qasm(head + "creg m[5];\nqreg q[2];\ncreg n[2];\nh q[1];");
        CHECK(classical.n == 2 && classical.size() == 1);

        // The registers together may use every qubit the IR allows.
        const Circuit full = parse_qasm(head + "qreg a[16];\nqreg b[16];\nx b[15];");
        validate(full);
        CHECK(full.n == 32 && full.gates[0].q[0] == 31);

        CHECK(parse_qasm("OPENQASM 2.0;\nqreg q[2];").size() == 0);
    }

    void test_qasm2_errors() {
        // Header and declarations.
        CHECK(fails_at("OPENQASM 2.0;\ninclude \"stdgates.inc\";", 2, 9));     // 3.0's library
        CHECK(fails_at("OPENQASM 2.0;\nqreg q[1];\nh q[0];", 3, 1));           // h before the include
        CHECK(fails_at("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n", 3, 1));     // no register
        CHECK(fails_at("OPENQASM 2.0;\nqreg q[0];", 2, 8));
        CHECK(fails_at("OPENQASM 2.0;\nqreg q[33];", 2, 8));
        CHECK(fails_at("OPENQASM 2.0;\nqreg a[20];\nqreg b[13];", 3, 8));      // 33 qubits in total
        CHECK(fails_at("OPENQASM 2.0;\nqreg q[2];\nqreg q[2];", 3, 6));        // same name twice
        CHECK(fails_at("OPENQASM 2.0;\nqreg [2] q;", 2, 6));
        CHECK(fails_at("OPENQASM 2.0;\nqreg q;", 2, 7));
        CHECK(fails_at("OPENQASM 2.0;\nqreg q[1];\ncreg c;", 3, 7));

        // Each version keeps its own declaration syntax.
        CHECK(fails_at("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqubit[2] q;", 3, 1));
        CHECK(fails_at("OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqreg q[2];", 3, 1));

        // Gates and statements outside 2.0 and qelib1.inc.
        CHECK(fails_at(kHeader2 + "mcx q[0],q[1],q[2];", 4, 1));
        CHECK(fails_at(kHeader2 + "csdg q[0],q[1];", 4, 1));
        CHECK(fails_at(kHeader2 + "swap q[0],q[1];", 4, 1));
        CHECK(fails_at(kHeader2 + "gphase(1);", 4, 1));
        CHECK(fails_at(kHeader2 + "measure q[0];", 4, 1));

        // Operands.
        CHECK(fails_at(kHeader2 + "x q[3];", 4, 5));
        CHECK(fails_at(kHeader2 + "x r[0];", 4, 3));
        CHECK(fails_at(kHeader2 + "cx q[0] q[1];", 4, 9));                     // missing ','
        const std::string two = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\nqreg flag[1];\n";
        CHECK(fails_at(two + "x flag[1];", 5, 8));

        CHECK(error_of(two + "x flag[1];") == "5:8: qubit index 1 is out of range for register 'flag' of size 1");
        CHECK(error_of("OPENQASM 2.0;\nqreg q[1];\nh q[0];") ==
              "3:1: gate 'h' is not defined; include \"qelib1.inc\" first");
        CHECK(error_of("OPENQASM 2.0;\nqreg a[20];\nqreg b[13];") ==
              "3:8: circuit would have more than 32 qubits");
        CHECK(error_of("OPENQASM 2.0;\nqreg q[2];\nqreg q[2];") == "3:6: register 'q' is already declared");
        CHECK(error_of(kHeader2 + "mcx q[0],q[1],q[2];") == "4:1: unknown gate or statement 'mcx'");
    }

    // ---- OpenQASM 2.0 gate definitions and qelib1.inc ----

    // Dense operators, column-major, as testgen::reference_unitary returns them.
    using Matrix = std::vector<cplx>;

    bool close(const Matrix& a, const Matrix& b) {
        if (a.size() != b.size()) return false;
        for (std::size_t i = 0; i < a.size(); ++i)
            if (std::abs(a[i] - b[i]) > 1e-12) return false;
        return true;
    }

    // The operator of a 2.0 program on qreg q[n].
    Matrix unitary2(int n, const std::string& body) {
        return testgen::reference_unitary(parse_qasm(
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[" + std::to_string(n) + "];\n" + body));
    }

    Matrix unitary_of(const Gate& g, int n) {
        Circuit c(n);
        c.add(g);
        return testgen::reference_unitary(c);
    }

    // V on q[1], controlled by q[0]: V acts when bit 0 of the basis index is set, that is on
    // indices 1 (target 0) and 3 (target 1).
    Matrix controlled(const Gate& v) {
        Matrix m(16, 0.0);
        const auto at = [&](int row, int col) -> cplx& { return m[col * 4 + row]; };
        at(0, 0) = 1;
        at(2, 2) = 1;
        at(1, 1) = v.m[0];
        at(1, 3) = v.m[1];
        at(3, 1) = v.m[2];
        at(3, 3) = v.m[3];
        return m;
    }

    // The 3-qubit identity with basis states i and j exchanged.
    Matrix swapped(std::size_t i, std::size_t j) {
        Matrix m(64, 0.0);
        for (std::size_t k = 0; k < 8; ++k) {
            const std::size_t to = k == i ? j : k == j ? i : k;
            m[k * 8 + to] = 1;
        }
        return m;
    }

    void test_qasm2_builtins() {
        const Circuit c = parse_qasm(kHeader2 + "U(0.3,0.2,0.1) q[0]; CX q[2],q[1];");
        CHECK(c.size() == 2);
        if (c.size() == 2) {
            CHECK(is_gate(c.gates[0], make_gate(K::U, 0, {0.3, 0.2, 0.1})));
            CHECK(is_gate(c.gates[1], make_gate(K::CX, 2, 1)));
        }
        CHECK(c.global_phase == 0.0);       // parser.hpp, point 3: 2.0's U is the IR's U

        // Built-ins need no include.
        const Circuit bare = parse_qasm("OPENQASM 2.0;\nqreg q[2];\nU(0,0,0) q[1]; CX q[0],q[1];");
        CHECK(bare.size() == 2);
    }

    // The qelib1.inc gates read from their definitions (parser.hpp, point 3).
    void test_qelib1_gates() {
        const Circuit u3 = parse_qasm(kHeader2 + "u3(0.3,0.2,0.1) q[1];");
        CHECK(u3.size() == 1 && is_gate(u3.gates[0], make_gate(K::U, 1, {0.3, 0.2, 0.1})));
        const Circuit u2 = parse_qasm(kHeader2 + "u2(0.2,0.1) q[1];");
        CHECK(u2.size() == 1 && is_gate(u2.gates[0], make_gate(K::U, 1, {kPi / 2, 0.2, 0.1})));
        const Circuit u1 = parse_qasm(kHeader2 + "u1(0.7) q[1];");
        CHECK(u1.size() == 1 && is_gate(u1.gates[0], make_gate(K::U, 1, {0, 0, 0.7})));

        CHECK(close(unitary2(1, "u1(0.7) q[0];"), unitary_of(make_gate(K::P, 0, {0.7, 0, 0}), 1)));
        CHECK(close(unitary2(1, "ry(0.7) q[0];"), unitary_of(make_gate(K::RY, 0, {0.7, 0, 0}), 1)));

        CHECK(close(unitary2(2, "cy q[0],q[1];"), controlled(make_gate(K::Y, 1))));
        CHECK(close(unitary2(2, "crz(0.7) q[0],q[1];"), controlled(make_gate(K::RZ, 1, {0.7, 0, 0}))));
        CHECK(close(unitary2(2, "cu1(0.7) q[0],q[1];"), controlled(make_gate(K::P, 1, {0.7, 0, 0}))));

        // parser.hpp, point 3: cu3 is controlled-U3; ch is controlled-H times e^{i pi/4}.
        CHECK(close(unitary2(2, "cu3(0.3,0.2,0.1) q[0],q[1];"),
                    controlled(make_gate(K::U, 1, {0.3, 0.2, 0.1}))));
        CHECK(close(unitary2(2, "cu3(2.9,-1.3,0.8) q[1],q[0];"),
                    unitary2(2, "cx q[0],q[1]; cx q[1],q[0]; cx q[0],q[1]; cu3(2.9,-1.3,0.8) q[0],q[1];"
                                " cx q[0],q[1]; cx q[1],q[0]; cx q[0],q[1];")));
        Matrix ch = controlled(make_gate(K::H, 1));
        for (cplx& z : ch) z *= std::polar(1.0, kPi / 4);
        CHECK(close(unitary2(2, "ch q[0],q[1];"), ch));

        // ccx a,b,c flips c when a and b are 1: 15 gates from the library's body.
        const Circuit ccx = parse_qasm(kHeader2 + "ccx q[0],q[1],q[2];");
        CHECK(ccx.size() == 15 && ccx.global_phase == 0.0);
        CHECK(close(unitary2(3, "ccx q[0],q[1],q[2];"), swapped(0b011, 0b111)));
        CHECK(close(unitary2(3, "ccx q[2],q[0],q[1];"), swapped(0b101, 0b111)));

        // The gates read directly keep their IR kinds.
        const Circuit direct = parse_qasm(kHeader2 + "h q[0]; cz q[0],q[1];");
        CHECK(direct.size() == 2 && direct.gates[0].kind == K::H && direct.gates[1].kind == K::CZ);
    }

    void test_gate_definitions() {
        // No parameters, then a call: the definition itself adds nothing.
        const Circuit bell = parse_qasm(kHeader2 + "gate bell a,b { h a; cx a,b; }\nbell q[1],q[2];");
        CHECK(bell.size() == 2);
        if (bell.size() == 2) {
            CHECK(is_gate(bell.gates[0], make_gate(K::H, 1)));
            CHECK(is_gate(bell.gates[1], make_gate(K::CX, 1, 2)));
        }

        // Parameters inside expressions, and the built-in U.
        const Circuit r = parse_qasm(kHeader2 + "gate r(t, u) a { rz(t/2) a; rx(-t*u) a; U(t, u, pi) a; }\n"
                                                "r(pi, 2) q[2];");
        CHECK(r.size() == 3);
        if (r.size() == 3) {
            CHECK(is_gate(r.gates[0], make_gate(K::RZ, 2, {kPi / 2, 0, 0})));
            CHECK(is_gate(r.gates[1], make_gate(K::RX, 2, {-kPi * 2, 0, 0})));
            CHECK(is_gate(r.gates[2], make_gate(K::U, 2, {kPi, 2, kPi})));
        }

        // A body may call earlier definitions, with its arguments in any order.
        const Circuit nested = parse_qasm(kHeader2 + "gate bell a,b { h a; cx a,b; }\n"
                                                     "gate rev a,b { bell b,a; }\nrev q[0],q[2];");
        CHECK(nested.size() == 2);
        if (nested.size() == 2) {
            CHECK(is_gate(nested.gates[0], make_gate(K::H, 2)));
            CHECK(is_gate(nested.gates[1], make_gate(K::CX, 2, 0)));
        }
        CHECK(close(unitary2(3, "gate tof a,b,c { ccx c,b,a; }\ntof q[0],q[1],q[2];"),
                    unitary2(3, "ccx q[2],q[1],q[0];")));

        // Empty parentheses, an empty body, and a call of each.
        CHECK(parse_qasm(kHeader2 + "gate g() a { x a; }\ng q[0]; g() q[1];").size() == 2);
        CHECK(parse_qasm(kHeader2 + "gate nop a { }\nnop q[0];").size() == 0);

        // A definition may come before the include, and may use a library name if the library
        // is never included.
        CHECK(parse_qasm("OPENQASM 2.0;\ngate g a { U(0,0,0) a; }\ninclude \"qelib1.inc\";\n"
                         "qreg q[1];\ng q[0]; h q[0];").size() == 2);
        const Circuit own_x = parse_qasm("OPENQASM 2.0;\nqreg q[1];\ngate x a { U(pi,0,pi) a; }\nx q[0];");
        CHECK(own_x.size() == 1 && own_x.gates[0].kind == K::U);

        // Comments and line breaks inside a definition.
        CHECK(parse_qasm(kHeader2 + "gate g(t) // a rotation\n  a\n{\n  rz(t) a; /* then */ h a;\n}\ng(1) q[0];")
                  .size() == 2);
    }

    void test_gate_definition_errors() {
        // A name is defined once: no redefinitions, of the file's own gates or of the library's.
        CHECK(fails_at(kHeader2 + "gate bell a,b { h a; }\ngate bell a { x a; }", 5, 6));
        CHECK(fails_at(kHeader2 + "gate ccx a,b,c { }", 4, 6));
        CHECK(fails_at(kHeader2 + "gate h a { }", 4, 6));
        CHECK(fails_at(kHeader2 + "gate U a { }", 4, 6));
        CHECK(fails_at(kHeader2 + "gate CX a,b { }", 4, 6));
        CHECK(fails_at("OPENQASM 2.0;\ngate ccx a,b,c { U(0,0,0) a; }\ninclude \"qelib1.inc\";", 3, 1));
        CHECK(fails_at("OPENQASM 2.0;\ngate h a { U(0,0,0) a; }\ninclude \"qelib1.inc\";", 3, 1));

        // The header.
        CHECK(fails_at(kHeader2 + "gate g(t,t) a { }", 4, 10));
        CHECK(fails_at(kHeader2 + "gate g a,a { }", 4, 10));
        CHECK(fails_at(kHeader2 + "gate g { }", 4, 8));
        CHECK(fails_at(kHeader2 + "gate g a h a;", 4, 10));
        CHECK(fails_at(kHeader2 + "gate g a { h a;", 4, 16));               // no closing '}'
        CHECK(fails_at(kHeader2 + "gate g(1) a { }", 4, 8));

        // The body: gates, parameters, and arguments.
        CHECK(fails_at(kHeader2 + "gate g a { foo a; }", 4, 12));
        CHECK(fails_at(kHeader2 + "gate g a { g a; }", 4, 12));             // not yet defined
        CHECK(fails_at(kHeader2 + "gate g a { h b; }", 4, 14));
        CHECK(fails_at(kHeader2 + "gate g a { h a[0]; }", 4, 15));
        CHECK(fails_at(kHeader2 + "gate g(t) a { rz(s) a; }", 4, 18));
        CHECK(fails_at(kHeader2 + "gate g a { cx a; }", 4, 12));
        CHECK(fails_at(kHeader2 + "gate g a,b { cx a,a; }", 4, 19));
        CHECK(fails_at(kHeader2 + "gate g a { rz a; }", 4, 12));
        CHECK(fails_at("OPENQASM 2.0;\ngate g a { h a; }", 2, 12));          // h needs the include

        // Calls.
        const std::string g = kHeader2 + "gate g(t) a,b { }\n";
        CHECK(fails_at(g + "g q[0],q[1];", 5, 1));
        CHECK(fails_at(g + "g(1) q[0];", 5, 1));
        CHECK(fails_at(g + "g(1) q[0],q[0];", 5, 11));
        CHECK(fails_at(kHeader2 + "ccx q[0],q[1],q[0];", 4, 15));
        CHECK(fails_at(kHeader2 + "gate g(t) a { rz(1/t) a; }\ng(0) q[0];", 5, 1));

        // OpenQASM 3.0 has no gate statements here.
        CHECK(fails_at(kHeader + "gate g a { }", 4, 1));

        CHECK(error_of(kHeader2 + "gate h a { }") == "4:6: gate 'h' is already defined");
        CHECK(error_of("OPENQASM 2.0;\ngate ccx a,b,c { U(0,0,0) a; }\ninclude \"qelib1.inc\";") ==
              "3:1: gate 'ccx' is already defined, so \"qelib1.inc\" cannot be included");
        CHECK(error_of(kHeader2 + "gate g a { h b; }") == "4:14: unknown qubit argument 'b'");
        CHECK(error_of(kHeader2 + "gate g(t) a { rz(s) a; }") == "4:18: unknown parameter 's'");
        CHECK(error_of(kHeader2 + "gate g a,b { cx a,a; }") == "4:19: gate 'cx' needs distinct qubits");
        CHECK(error_of(kHeader2 + "gate g(t) a { rz(1/t) a; }\ng(0) q[0];") ==
              "5:1: gate 'g': a parameter of 'g' does not evaluate to a finite number");
    }

    void test_messages() {
        CHECK(error_of("OPENQASM 4.0;") == "1:10: unsupported version '4.0'; only 2.0 and 3.0 are supported");
        CHECK(error_of(kHeader + "cx q[0] q[1];") == "4:9: expected ';' but found 'q'");
        CHECK(error_of(kHeader + "h q[0]") == "4:7: expected ';' but found end of input");
        CHECK(error_of(kHeader + "h q[3];") == "4:5: qubit index 3 is out of range for register 'q' of size 3");
        CHECK(error_of(kHeader + "rz(1,2) q[0];") == "4:1: gate 'rz' takes 1 parameter but got 2");
        CHECK(error_of(kHeader + "h q[0],q[1];") == "4:1: gate 'h' acts on 1 qubit but got 2");
        CHECK(error_of(kHeader + "cx q[0];") == "4:1: gate 'cx' acts on 2 qubits but got 1");
        CHECK(error_of(kHeader + "foo q[0];") == "4:1: unknown gate or statement 'foo'");
        CHECK(error_of("OPENQASM 3.0;\ninclude \"qelib1.inc\";") ==
              "2:9: cannot include \"qelib1.inc\"; only \"stdgates.inc\" is supported");
    }
} // namespace

int main() {
    test_program();
    test_every_gate();
    test_global_phase();
    test_numbers();
    test_expressions();
    test_expression_errors();
    test_layout();
    test_comments();
    test_dag();
    test_header_errors();
    test_gate_errors();
    test_qasm2_program();
    test_qasm2_gates();
    test_qasm2_registers();
    test_qasm2_errors();
    test_qasm2_builtins();
    test_qelib1_gates();
    test_gate_definitions();
    test_gate_definition_errors();
    test_messages();
    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("all parser tests passed\n");
    return 0;
}
