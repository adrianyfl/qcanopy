#include "qcanopy/dataset_io.hpp"
#include <charconv>
#include <cmath>
#include <iostream>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>

using namespace qcanopy::testgen;
namespace {
const char* help = R"(qcanopy_generate --format qasm|ir --output DIR [options]
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
)";
uint64_t integer(const std::string& s) {
    uint64_t out = 0;
    const auto r = std::from_chars(s.data(),s.data()+s.size(),out);
    if (s.empty() || r.ec != std::errc{} || r.ptr != s.data()+s.size()) throw std::invalid_argument("expected decimal UINT64: " + s);
    return out;
}
std::size_t size(const std::string& s) {
    const auto n = integer(s);
    if (n > std::numeric_limits<std::size_t>::max()) throw std::invalid_argument("integer exceeds platform size_t");
    return static_cast<std::size_t>(n);
}
double real(const std::string& s) {
    std::istringstream in(s); in.imbue(std::locale::classic());
    double out; in >> std::noskipws >> out;
    if (!in || !in.eof() || !std::isfinite(out)) throw std::invalid_argument("expected finite number: " + s);
    return out;
}
}
int main(int argc, char** argv) {
    try {
        GeneratorConfig c; OutputOptions o; std::set<std::string> seen;
        for (int i = 1; i < argc; ++i) {
            const std::string key = argv[i];
            if (key == "--help") { std::cout << help; return 0; }
            if (!seen.insert(key).second) throw std::invalid_argument("duplicate option: " + key);
            if (i+1 == argc) throw std::invalid_argument("missing value for " + key);
            const std::string v = argv[++i];
            if (key == "--format") o.format = v;
            else if (key == "--output") o.destination = v;
            else if (key == "--qubits") {
                const auto n = integer(v);
                if (n > qcanopy::kMaxQubits) throw std::invalid_argument("qubits exceeds 32");
                c.qubits = static_cast<int>(n);
            }
            else if (key == "--gates") c.gates = size(v);
            else if (key == "--pairs") c.pairs = size(v);
            else if (key == "--prefix-group-size") c.prefix_group_size = size(v);
            else if (key == "--seed") c.seed = integer(v);
            else if (key == "--prefix-fraction") c.prefix_fraction = real(v);
            else if (key == "--equivalent-fraction") c.equivalent_fraction = real(v);
            else if (key == "--gate-set") c.gate_set = v;
            else if (key == "--angle-mode") c.angle_mode = v;
            else if (key == "--verify") c.verify = v;
            else throw std::invalid_argument("unknown option: " + key);
        }
        validate_config(c); validate_output_options(o);
        const auto cases = generate_batch(c);
        write_dataset(o,c,cases);
        std::size_t eq = 0, verified = 0;
        for (const auto& p : cases) { eq += p.expected == ExpectedRelation::Equivalent; verified += p.metadata.verification == "passed"; }
        std::cout << "Generated " << cases.size() << " pairs (" << eq << " equivalent, " << cases.size()-eq
                  << " non-equivalent), " << 2*cases.size() << " " << o.format << " circuits; " << verified
                  << " numerically verified.\nManifest: " << (o.destination/"manifest.jsonl").string() << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "qcanopy_generate: " << e.what() << "\nUse --help for usage.\n"; return 1;
    }
}
