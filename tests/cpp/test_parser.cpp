// Tests for parse_qasm (qcanopy/parser.hpp).
#include "qcanopy/parser.hpp"

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

    void test_numbers() {
        const auto angle = [](const std::string& text) {
            const Circuit c = parse_body("rz(" + text + ") q[0];");
            return c.size() == 1 ? c.gates[0].params[0] : -999.0;
        };
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
        CHECK(fails_at("OPENQASM 2.0;\nqreg q[2];", 1, 10));
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
        CHECK(body_fails_at("rz(--1) q[0];", 4, 5));
        CHECK(body_fails_at("rz(q) q[0];", 4, 4));
        CHECK(body_fails_at("rz(pi) q[0];", 4, 4));                         // no named constants
        CHECK(body_fails_at("rz(1/2) q[0];", 4, 5));                        // no arithmetic
        CHECK(body_fails_at("rz(+1) q[0];", 4, 4));
        CHECK(body_fails_at("rz(1e999) q[0];", 4, 4));

        CHECK(body_fails_at("gphase 1;", 4, 8));
        CHECK(body_fails_at("gphase(1) q[0];", 4, 11));
        CHECK(body_fails_at("gphase();", 4, 8));

        // Lexer errors come through unchanged.
        CHECK(body_fails_at("h q[0]; @", 4, 9));
        CHECK(body_fails_at("rz(1e) q[0];", 4, 4));
    }

    void test_messages() {
        CHECK(error_of("OPENQASM 2.0;") == "1:10: unsupported version '2.0'; only 3.0 is supported");
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
    test_layout();
    test_comments();
    test_dag();
    test_header_errors();
    test_gate_errors();
    test_messages();
    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("all parser tests passed\n");
    return 0;
}
