// Tests for `qcanopy check`: argument parsing, stimulus building, and check_circuits().
#include "check.hpp"

#include "qcanopy/parser.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
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

    void test_proved_by_rewrite() {
        // S S fuses to Z during normalize; the rest peels off.
        const Circuit a = qasm(2, "h q[0];\ns q[1];\ns q[1];\ncx q[0],q[1];\n");
        const Circuit b = qasm(2, "h q[0];\nz q[1];\ncx q[0],q[1];\n");
        const CheckOutcome r = check_circuits(a, b, CheckOptions{});
        CHECK(r.equiv.equivalent && r.proved_by_rewrite);
        CHECK(r.residual_a == 0 && r.residual_b == 0 && r.equiv.stimuli_checked == 0);
        CHECK(std::abs(phase_of(r.equiv)) < 1e-9);
    }

    void test_global_phase_from_rewrite() {
        // gphase is carried through normalize + peel into the reported phase.
        const Circuit a = qasm(1, "gphase(0.25);\nh q[0];\n");
        const Circuit b = qasm(1, "h q[0];\n");
        const CheckOutcome r = check_circuits(a, b, CheckOptions{});
        CHECK(r.equiv.equivalent && r.proved_by_rewrite);
        CHECK(std::abs(phase_of(r.equiv) - 0.25) < 1e-9);
    }

    void test_needs_simulation() {
        // CX with control and target swapped by Hadamard conjugation: nothing peels.
        const Circuit a = qasm(2, "cx q[0],q[1];\n");
        const Circuit b = qasm(2, "h q[0];\nh q[1];\ncx q[1],q[0];\nh q[0];\nh q[1];\n");
        for (const bool rewrite : {true, false}) {
            CheckOptions o;
            o.rewrite = rewrite;
            const CheckOutcome r = check_circuits(a, b, o);
            CHECK(r.equiv.equivalent && !r.proved_by_rewrite);
            CHECK(r.stimuli == 16 && r.equiv.stimuli_checked == 16);
        }
    }

    void test_not_equivalent() {
        const Circuit a = qasm(2, "h q[0];\ncx q[0],q[1];\nz q[0];\n");
        const Circuit b = qasm(2, "h q[0];\ncx q[0],q[1];\n");
        CheckOptions o;
        o.stimuli = {"all-basis"};
        const CheckOutcome r = check_circuits(a, b, o);
        CHECK(!r.equiv.equivalent && !r.proved_by_rewrite);
        CHECK(r.residual_a == 1 && r.residual_b == 0);  // only the Z is left to simulate
        CHECK(r.equiv.mismatch != Mismatch::None);
    }

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
    test_build_stimulus();
    test_proved_by_rewrite();
    test_global_phase_from_rewrite();
    test_needs_simulation();
    test_not_equivalent();
    test_width_mismatch();
    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("all tests passed\n");
    return 0;
}
