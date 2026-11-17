// Round-trip tests for the parser: circuits from the generator are written with write_qasm
// and read back with parse_qasm (qcanopy/parser.hpp).
#include "qcanopy/dataset_io.hpp"
#include "qcanopy/parser.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace qcanopy;
using namespace qcanopy::testgen;

namespace {
    int failures = 0;

#define CHECK(cond) \
    do { if (!(cond)) { std::printf("%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)

    using K = GateKind;

    constexpr double kPi = 3.14159265358979323846;

    // The writer subtracts theta/2 per U gate and the parser adds it back, in a different
    // order, so the global phase returns only up to rounding. Gates return exactly.
    constexpr double kPhaseTol = 1e-10;

    std::string to_qasm(const Circuit& c) {
        std::ostringstream out;
        write_qasm(out, c);
        return out.str();
    }

    // What write_qasm's output reads back as. The writer has no sx/sxdg: it emits
    // rx(+-pi/2) and moves the remaining +-pi/4 into gphase. Every other gate is unchanged.
    Circuit expected(const Circuit& c) {
        Circuit e(c.n);
        e.global_phase = c.global_phase;
        for (const Gate& g : c.gates) {
            if (g.kind == K::SX) {
                e.add(make_gate(K::RX, g.q[0], {kPi / 2, 0, 0}));
                e.global_phase += kPi / 4;
            } else if (g.kind == K::SXdg) {
                e.add(make_gate(K::RX, g.q[0], {-kPi / 2, 0, 0}));
                e.global_phase -= kPi / 4;
            } else {
                e.add(g);
            }
        }
        return e;
    }

    bool same_circuit(const Circuit& got, const Circuit& want) {
        if (got.n != want.n || got.size() != want.size()) return false;
        if (std::abs(got.global_phase - want.global_phase) > kPhaseTol) return false;
        for (std::size_t i = 0; i < got.size(); ++i) {
            const Gate& a = got.gates[i];
            const Gate& b = want.gates[i];
            if (!same_gate(a, b, 0.0) || a.q != b.q || a.params != b.params) return false;
        }
        return true;
    }

    // Full operators agree entry by entry, global phase included. Needs n <= kDenseQubitCap.
    bool same_unitary(const Circuit& a, const Circuit& b) {
        const std::vector<cplx> ua = reference_unitary(a), ub = reference_unitary(b);
        if (ua.size() != ub.size()) return false;
        for (std::size_t i = 0; i < ua.size(); ++i)
            if (std::abs(ua[i] - ub[i]) > 1e-9) return false;
        return true;
    }

    // Writes and re-reads both circuits of every pair in a generated batch.
    void roundtrip_batch(const GeneratorConfig& config) {
        std::size_t bad = 0;
        for (const CircuitPair& pair : generate_batch(config)) {
            for (const Circuit* c : {&pair.a, &pair.b}) {
                const Circuit got = parse_qasm(to_qasm(*c));
                validate_circuit(got);
                bool ok = same_circuit(got, expected(*c));
                if (ok && c->n <= kDenseQubitCap) ok = same_unitary(got, *c);
                if (!ok) ++bad;
            }
        }
        if (bad)
            std::printf("round trip failed for %zu circuit(s): qubits=%d gates=%zu gate_set=%s angle_mode=%s seed=%llu\n",
                        bad, config.qubits, config.gates, config.gate_set.c_str(), config.angle_mode.c_str(),
                        static_cast<unsigned long long>(config.seed));
        CHECK(bad == 0);
    }

    void test_generated_batches() {
        roundtrip_batch(GeneratorConfig{});                     // the generator's defaults

        for (const char* gate_set : {"clifford", "clifford-t", "parameterized", "mixed"})
            for (const char* angle_mode : {"discrete", "continuous"}) {
                GeneratorConfig c;
                c.pairs = 8;
                c.seed = 11;
                c.gate_set = gate_set;
                c.angle_mode = angle_mode;
                roundtrip_batch(c);
            }

        // Circuit widths from the narrowest to the widest the IR allows. Above
        // kDenseQubitCap only the gate lists are compared.
        for (int qubits : {1, 2, 3, 5, 8, 16, 32}) {
            GeneratorConfig c;
            c.qubits = qubits;
            c.gates = 300;
            c.pairs = 4;
            c.seed = 100 + static_cast<uint64_t>(qubits);
            roundtrip_batch(c);
        }
    }

    // One circuit holding every named gate kind, including the ones the generator's default
    // gate sets may not draw.
    void test_every_kind() {
        Circuit c(2);
        c.global_phase = 0.3;
        for (int i = 0; i < static_cast<int>(K::U1Q); ++i) {
            const K k = static_cast<K>(i);
            std::array<double, 3> p{};
            for (int j = 0; j < info(k).nparams; ++j) p[j] = 0.25 * (j + 1);
            c.add(info(k).nq == 1 ? make_gate(k, 0, p) : make_gate(k, 1, 0));
        }
        const Circuit got = parse_qasm(to_qasm(c));
        CHECK(got.size() == c.size());
        CHECK(same_circuit(got, expected(c)));
        CHECK(same_unitary(got, c));
    }

    void test_dag() {
        GeneratorConfig config;
        config.pairs = 4;
        for (const CircuitPair& pair : generate_batch(config)) {
            const std::string text = to_qasm(pair.a);
            const CircuitDAG dag = parse_qasm_dag(text);
            const Circuit want = expected(pair.a);
            CHECK(dag.n == want.n);
            CHECK(dag.nodes.size() == want.size());
            // Walking the DAG gives the same circuit as building the DAG by hand.
            CHECK(same_circuit(dag.to_circuit(), CircuitDAG(want).to_circuit()));
            CHECK(same_unitary(dag.to_circuit(), pair.a));
        }
    }

    void test_files() {
        namespace fs = std::filesystem;
        const fs::path root = fs::temp_directory_path() /
            ("qcanopy-parser-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        fs::create_directories(root);

        const Circuit c = generate_circuit(GeneratorConfig{}, 42);
        const std::string good = (root / "good.qasm").string();
        { std::ofstream out(good, std::ios::binary); out << to_qasm(c); }

        CHECK(same_circuit(parse_qasm_file(good), expected(c)));
        const CircuitDAG dag = parse_qasm_file_dag(good);
        CHECK(dag.n == c.n && dag.nodes.size() == c.size());

        // A parse error names the file and keeps its position.
        const std::string bad = (root / "bad.qasm").string();
        { std::ofstream out(bad, std::ios::binary); out << "OPENQASM 3.0;\nqubit[2] q;\nh q[0];\n"; }
        try {
            (void)parse_qasm_file(bad);
            CHECK(false);
        } catch (const ParseError& e) {
            CHECK(e.line() == 3 && e.column() == 1);
            CHECK(std::string(e.what()) ==
                  bad + ":3:1: gate 'h' is not defined; include \"stdgates.inc\" first");
        }
        try {
            (void)parse_qasm_file_dag(bad);
            CHECK(false);
        } catch (const ParseError&) {}

        // A missing file is an I/O failure, not a parse error.
        const std::string missing = (root / "missing.qasm").string();
        try {
            (void)parse_qasm_file(missing);
            CHECK(false);
        } catch (const ParseError&) {
            CHECK(false);
        } catch (const std::runtime_error& e) {
            CHECK(std::string(e.what()) == "cannot open " + missing);
        }

        fs::remove_all(root);
    }
} // namespace

int main() {
    test_generated_batches();
    test_every_kind();
    test_dag();
    test_files();
    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("all parser round-trip tests passed\n");
    return 0;
}
