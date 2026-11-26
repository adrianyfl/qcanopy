// `qcanopy generate`: labeled circuit-pair datasets (generator.hpp, dataset_io.hpp).
// The older qcanopy_generate tool runs this too; options and output are unchanged.

#include "commands.hpp"
#include "common.hpp"

#include "qcanopy/dataset_io.hpp"

#include <iostream>
#include <limits>
#include <string>

namespace qcanopy::cli {

    namespace {
        using namespace qcanopy::testgen;

        const char* const help = R"(qcanopy generate --format qasm|ir --output DIR [options]
Generate labeled circuit pairs with planted prefixes and normalized sharing metrics.
  --format qasm|ir          Required. OpenQASM 3.0 + stdgates.inc or JSON IR v1.
  --output DIR              Required. Absent or empty destination; never overwrite.
  --qubits N                1..32 (default 4)
  --gates L                 Base gate count, not depth (default 100; zero allowed)
  --pairs B                 Positive pair count (default 32)
  --prefix-fraction R       Planted floor(R*L) common gates, R in [0,1] (default .5)
  --equivalent-fraction E   floor(E*B+.5) equivalents, seeded shuffle (default .5)
  --seed UINT64             Decimal unsigned master seed (default 759)
  --gate-set SET            clifford|clifford-t|parameterized|mixed (default mixed)
  --angle-mode MODE         discrete|continuous|mixed (default mixed)
  --prefix-group-size G     Positive consecutive sharing group size (default 1)
  --verify small|none       Full-operator + normalization checks up to 5 qubits
                            (default small); larger cases retain construction proofs.
  --help                    Show this help.
Equivalent rewrites can change gate counts and measured sharing. Outputs retain
raw circuits, global phase, labels, seeds, proofs, metrics, and manifest.jsonl.
Check a generated pair with: qcanopy check DIR/cases/case-000000.a.qasm DIR/cases/case-000000.b.qasm
)";

        std::size_t parse_size(const std::string& key, const std::string& s) {
            const uint64_t n = parse_uint(key, s);
            if (n > std::numeric_limits<std::size_t>::max()) throw UsageError(key + " exceeds platform size_t");
            return static_cast<std::size_t>(n);
        }
    } // namespace

    int generate_main(int argc, const char* const* argv) {
        return guarded("generate", [&] {
            const Args args(argc, argv, {},
                            {"--format", "--output", "--qubits", "--gates", "--pairs", "--prefix-group-size", "--seed",
                             "--prefix-fraction", "--equivalent-fraction", "--gate-set", "--angle-mode", "--verify"});
            if (args.help()) { std::cout << help; return kExitOk; }
            if (!args.positional().empty()) throw UsageError("unexpected argument '" + args.positional()[0] + "'");

            GeneratorConfig c;
            OutputOptions o;
            if (auto v = args.get("--format")) o.format = *v;
            if (auto v = args.get("--output")) o.destination = *v;
            if (auto v = args.get("--qubits")) {
                const uint64_t n = parse_uint("--qubits", *v);
                if (n > static_cast<uint64_t>(kMaxQubits)) throw UsageError("--qubits exceeds " + std::to_string(kMaxQubits));
                c.qubits = static_cast<int>(n);
            }
            if (auto v = args.get("--gates")) c.gates = parse_size("--gates", *v);
            if (auto v = args.get("--pairs")) c.pairs = parse_size("--pairs", *v);
            if (auto v = args.get("--prefix-group-size")) c.prefix_group_size = parse_size("--prefix-group-size", *v);
            if (auto v = args.get("--seed")) c.seed = parse_uint("--seed", *v);
            if (auto v = args.get("--prefix-fraction")) c.prefix_fraction = parse_real("--prefix-fraction", *v);
            if (auto v = args.get("--equivalent-fraction")) c.equivalent_fraction = parse_real("--equivalent-fraction", *v);
            if (auto v = args.get("--gate-set")) c.gate_set = *v;
            if (auto v = args.get("--angle-mode")) c.angle_mode = *v;
            if (auto v = args.get("--verify")) c.verify = *v;

            validate_config(c);
            validate_output_options(o);
            const auto cases = generate_batch(c);
            write_dataset(o, c, cases);

            std::size_t eq = 0, verified = 0;
            for (const auto& p : cases) {
                eq += p.expected == ExpectedRelation::Equivalent;
                verified += p.metadata.verification == "passed";
            }
            std::cout << "Generated " << cases.size() << " pairs (" << eq << " equivalent, " << cases.size() - eq
                      << " non-equivalent), " << 2 * cases.size() << " " << o.format << " circuits; " << verified
                      << " numerically verified.\nManifest: " << (o.destination / "manifest.jsonl").string() << '\n';
            return kExitOk;
        });
    }

} // namespace qcanopy::cli
