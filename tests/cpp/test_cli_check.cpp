// Tests for `qcanopy check`: argument parsing, stimulus building, and check_circuits().
#include "check.hpp"

#include "qcanopy/parser.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

using namespace qcanopy;
using namespace qcanopy::cli;

namespace {
    int failures = 0;

#define CHECK(cond) \
    do { if (!(cond)) { std::printf("%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)

    CheckOptions parse(std::vector<const char*> args) {
        args.insert(args.begin(), "check");
        return parse_check_args(static_cast<int>(args.size()), args.data());
    }

    bool rejects(std::vector<const char*> args) {
        try { parse(std::move(args)); } catch (const UsageError&) { return true; }
        return false;
    }

    Circuit qasm(int n, const std::string& body) {
        return parse_qasm("OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[" + std::to_string(n) + "] q;\n" + body);
    }

    double phase_of(const EquivResult& e) { return std::arg(e.global_phase); }

    void test_parse_defaults() {
        const CheckOptions o = parse({"a.qasm", "b.qasm"});
        CHECK(o.file_a == "a.qasm" && o.file_b == "b.qasm");
        CHECK((o.stimuli == std::vector<std::string>{"haar", "product"}));
        CHECK(o.count == 8 && o.tol == 1e-9 && o.rewrite && !o.json && !o.quiet && !o.help);
    }

    void test_parse_options() {
        const CheckOptions o = parse({"--stimulus", "all-basis,haar,haar", "a", "--count", "3", "b",
                                      "--seed", "42", "--tol", "1e-6", "--no-rewrite", "--json"});
        CHECK(o.file_a == "a" && o.file_b == "b");
        CHECK((o.stimuli == std::vector<std::string>{"all-basis", "haar"}));  // duplicates dropped
        CHECK(o.count == 3 && o.seed == 42 && o.tol == 1e-6 && !o.rewrite && o.json);
        CHECK(parse({"--help"}).help);
        CHECK(parse({"a", "--help", "--bogus"}).help);
        CHECK(parse({"a", "b", "--overlap"}).rewriter.overlap);
    }

    void test_parse_errors() {
        CHECK(rejects({}));
        CHECK(rejects({"a"}));
        CHECK(rejects({"a", "b", "c"}));
        CHECK(rejects({"a", "b", "--bogus"}));
        CHECK(rejects({"a", "b", "--count"}));
        CHECK(rejects({"a", "b", "--count", "0"}));
        CHECK(rejects({"a", "b", "--count", "-1"}));
        CHECK(rejects({"a", "b", "--seed", "12x"}));
        CHECK(rejects({"a", "b", "--tol", "0"}));
        CHECK(rejects({"a", "b", "--tol", "nan"}));
        CHECK(rejects({"a", "b", "--stimulus", "haar,"}));
        CHECK(rejects({"a", "b", "--stimulus", "gaussian"}));
        CHECK(rejects({"a", "b", "--json", "--quiet"}));
        CHECK(rejects({"a", "b", "--json", "--json"}));
        CHECK(rejects({"a", "b", "--overlap", "--no-rewrite"}));
    }

    // ---- Args: the option parser every subcommand uses ----

    Args args(std::vector<const char*> v) {
        v.insert(v.begin(), "cmd");
        return Args(static_cast<int>(v.size()), v.data(), {"--flag"}, {"--key"});
    }

    bool args_reject(std::vector<const char*> v) {
        try { args(std::move(v)); } catch (const UsageError&) { return true; }
        return false;
    }

    void test_args() {
        const Args a = args({"x", "--key", "-1", "--flag", "y"});
        CHECK((a.positional() == std::vector<std::string>{"x", "y"}));
        CHECK(a.has("--flag") && a.get("--key") == std::string("-1") && !a.get("--other"));
        CHECK(!a.help() && args({"x", "-h", "--bogus"}).help());
        CHECK(args_reject({"--bogus"}));
        CHECK(args_reject({"--flag", "--flag"}));
        CHECK(args_reject({"--key", "1", "--key", "2"}));
        CHECK(args_reject({"--key"}));
        bool threw = false;
        try { a.expect_positional(3, 3, "three files"); } catch (const UsageError&) { threw = true; }
        CHECK(threw);
    }

    void test_build_stimulus() {
        CheckOptions o;
        o.stimuli = {"haar", "basis", "all-basis"};
        o.count = 3;
        CHECK(build_stimulus(2, o).size() == 3 + 3 + 4);
        o.stimuli = {"all-basis"};
        bool threw = false;
        try { build_stimulus(kMaxAllBasisQubits + 1, o); } catch (const UsageError&) { threw = true; }
        CHECK(threw);
    }

    // A config whose interpreter does not exist, so the rewriter always fails.
    CheckOptions without_rewriter() {
        CheckOptions o;
        o.rewriter.python = "/nonexistent/qcanopy-test-python";
        return o;
    }

    // ---- parse_pair_proof: reading the prover's output ----

    const char* const kInconclusive = R"({"schema":"qrewrite-pair/1","swap_removal":{},
        "layout":{"cores_used":true,"policy":"x"},
        "proof":{"verdict":"inconclusive","proof":[{"kind":"structural"}],
          "residual_left":{"schema":"qdag/1","num_qubits":2,"global_phase":0.5,
            "nodes":[{"id":"g9","op":"cx","qubits":[0,1],"params":[]},
                     {"id":"g2","op":"rz","qubits":[0],"params":[0.25]}],
            "edges":[{"source":"g2","target":"g9","qubit":0}]},
          "residual_right":{"schema":"qdag/1","num_qubits":2,"global_phase":0,"nodes":[],"edges":[]}}})";

    std::string with(std::string text, const std::string& from, const std::string& to) {
        text.replace(text.find(from), from.size(), to);
        return text;
    }

    bool rejects_proof(const std::string& text) {
        std::istringstream in(text);
        try { parse_pair_proof(in, 2); } catch (const RewriterError&) { return true; }
        return false;
    }

    void test_parse_pair_proof() {
        std::istringstream in(kInconclusive);
        const PairProof p = parse_pair_proof(in, 2);
        CHECK(!p.proved && p.cores_used && p.steps == 1);
        // Nodes are listed out of order; the edge puts rz before cx.
        CHECK(p.left.size() == 2 && p.left.gates[0].kind == GateKind::RZ && p.left.gates[1].kind == GateKind::CX);
        CHECK(p.left.gates[0].params[0] == 0.25 && p.left.global_phase == 0.0);
        CHECK(p.right.size() == 0);

        CHECK(rejects_proof(with(kInconclusive, "\"rz\"", "\"rzz\"")));        // no IR gate
        CHECK(rejects_proof(with(kInconclusive, "[0.25]", "[\"pi/2\"]")));      // symbolic parameter
        CHECK(rejects_proof(with(kInconclusive, "\"num_qubits\":2", "\"num_qubits\":3")));
        CHECK(rejects_proof(with(kInconclusive, "\"inconclusive\"", "\"maybe\"")));
        CHECK(rejects_proof(with(kInconclusive, "\"inconclusive\"", "\"equivalent_up_to_global_phase\"")));
        CHECK(rejects_proof("{not json"));
    }

    // ---- check_circuits without the rewriter (always available) ----

    void test_fallback_when_rewriter_fails() {
        const Circuit a = qasm(2, "cx q[0],q[1];\n");
        const Circuit b = qasm(2, "h q[0];\nh q[1];\ncx q[1],q[0];\nh q[0];\nh q[1];\n");
        const CheckOutcome r = check_circuits(a, b, without_rewriter());
        CHECK(r.rewrite == RewriteStatus::Failed && !r.rewrite_error.empty());
        CHECK(r.equiv.equivalent && !r.phase_measured);
        CHECK(r.residual_a == 1 && r.residual_b == 5 && r.equiv.stimuli_checked == 16);
    }

    void test_no_rewrite() {
        CheckOptions o;
        o.rewrite = false;
        o.stimuli = {"all-basis"};
        const CheckOutcome r = check_circuits(qasm(1, "gphase(0.25);\nh q[0];\n"), qasm(1, "h q[0];\n"), o);
        CHECK(r.rewrite == RewriteStatus::Skipped && r.equiv.equivalent && r.stimuli == 2);
        CHECK(std::abs(phase_of(r.equiv) - 0.25) < 1e-9);
    }

    // ---- check_circuits with the real prover (needs Python 3.9+ at configure time) ----

#if QCANOPY_TEST_HAVE_PYTHON
    void test_rewriter_proves() {
        // rz(0.3) rz(-0.3) is an inverse pair; the rest matches structurally.
        const Circuit a = qasm(2, "gphase(0.25);\nh q[0];\nrz(0.3) q[1];\nrz(-0.3) q[1];\ncx q[0],q[1];\n");
        const Circuit b = qasm(2, "h q[0];\ncx q[0],q[1];\n");
        const CheckOutcome r = check_circuits(a, b, CheckOptions{});
        CHECK(r.rewrite == RewriteStatus::Proved && r.proof_steps > 0);
        CHECK(r.equiv.equivalent && r.residual_a == 0 && r.residual_b == 0 && r.stimuli == 0);
        CHECK(r.phase_measured && std::abs(phase_of(r.equiv) - 0.25) < 1e-9);  // phase from the full run
    }

    void test_rewriter_residual_not_equivalent() {
        const Circuit a = qasm(2, "h q[0];\ncx q[0],q[1];\nt q[0];\n");
        const Circuit b = qasm(2, "h q[0];\ncx q[0],q[1];\nz q[0];\n");
        CheckOptions o;
        o.stimuli = {"all-basis"};
        const CheckOutcome r = check_circuits(a, b, o);
        CHECK(r.rewrite == RewriteStatus::Inconclusive);
        CHECK(r.residual_a == 1 && r.residual_b == 1);  // only t vs z is simulated
        CHECK(!r.equiv.equivalent && !r.rewriter_contradicted && !r.phase_measured);
    }

    void test_rewriter_residual_equivalent() {
        // s s vs z: the prover may or may not prove this, but the verdict and phase must hold.
        const Circuit a = qasm(2, "gphase(-0.5);\nh q[0];\ns q[1];\ns q[1];\ncx q[0],q[1];\n");
        const Circuit b = qasm(2, "h q[0];\nz q[1];\ncx q[0],q[1];\n");
        const CheckOutcome r = check_circuits(a, b, CheckOptions{});
        CHECK(r.rewrite == RewriteStatus::Proved || r.rewrite == RewriteStatus::Inconclusive);
        CHECK(r.residual_a < a.size());
        CHECK(r.equiv.equivalent && r.phase_measured && !r.rewriter_contradicted);
        CHECK(std::abs(phase_of(r.equiv) + 0.5) < 1e-9);
    }

    void test_rewriter_overlap() {
        CheckOptions o;
        o.rewriter.overlap = true;
        const CheckOutcome r = check_circuits(qasm(2, "cx q[0],q[1];\n"), qasm(2, "cz q[0],q[1];\n"), o);
        CHECK(r.rewrite != RewriteStatus::Skipped && !r.equiv.equivalent);
    }
#endif

    void test_width_mismatch() {
        bool threw = false;
        try { check_circuits(qasm(1, ""), qasm(2, ""), CheckOptions{}); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }
}

int main() {
    test_parse_defaults();
    test_parse_options();
    test_parse_errors();
    test_args();
    test_build_stimulus();
    test_parse_pair_proof();
    test_fallback_when_rewriter_fails();
    test_no_rewrite();
#if QCANOPY_TEST_HAVE_PYTHON
    test_rewriter_proves();
    test_rewriter_residual_not_equivalent();
    test_rewriter_residual_equivalent();
    test_rewriter_overlap();
#else
    std::printf("Python 3.9+ not found at configure time: skipping rewriter tests\n");
#endif
    test_width_mismatch();
    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("all tests passed\n");
    return 0;
}
