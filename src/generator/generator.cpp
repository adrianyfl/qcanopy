#include "qcanopy/generator.hpp"
#include "qcanopy/dag.hpp"
#include "qcanopy/normalize.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace qcanopy::testgen {
namespace {
constexpr double pi = 3.14159265358979323846;
// SplitMix64, with unsigned wraparound defined by C++.
uint64_t mix(uint64_t z) {
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}
uint64_t derive(uint64_t master, uint64_t domain, uint64_t index) {
    return mix(master ^ mix(domain) ^ mix(index + UINT64_C(0x9e3779b97f4a7c15)));
}
struct Random {
    uint64_t state;
    uint64_t next() { state += UINT64_C(0x9e3779b97f4a7c15); return mix(state); }
    uint64_t bounded(uint64_t n) {
        const uint64_t threshold = (uint64_t{0} - n) % n;
        uint64_t x;
        do { x = next(); } while (x < threshold);
        return x % n;
    }
    double unit() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
};
double angle(Random& r, const std::string& mode) {
    bool discrete = mode == "discrete";
    if (mode == "mixed") discrete = r.bounded(2) == 0;
    if (discrete) return (static_cast<int>(r.bounded(17)) - 8) * (pi / 4);
    // A bounded offset from pi/4 grid points avoids normalization boundary noise.
    double x = (2 * r.unit() - 1) * (2 * pi);
    const double nearest = std::round(x / (pi / 4)) * (pi / 4);
    if (std::abs(x - nearest) < 1e-6) x += 1e-5;
    return x;
}
std::vector<GateKind> pool(const GeneratorConfig& c) {
    using K = GateKind;
    std::vector<K> kinds;
    if (c.gate_set != "parameterized") kinds = {K::H,K::X,K::Y,K::Z,K::S,K::Sdg,K::SX,K::SXdg};
    if (c.gate_set == "clifford-t" || c.gate_set == "mixed") {
        kinds.push_back(K::T); kinds.push_back(K::Tdg);
    }
    if (c.gate_set == "parameterized" || c.gate_set == "mixed")
        for (auto k : {K::RX,K::RY,K::RZ,K::P,K::U}) kinds.push_back(k);
    if (c.qubits > 1) for (auto k : {K::CX,K::CZ,K::SWAP}) kinds.push_back(k);
    return kinds;
}
Circuit random_circuit(const GeneratorConfig& c, std::size_t count, uint64_t seed) {
    Random r{seed};
    Circuit out(c.qubits);
    const auto kinds = pool(c);
    out.gates.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const auto kind = kinds[static_cast<std::size_t>(r.bounded(kinds.size()))];
        const auto q = static_cast<uint8_t>(r.bounded(static_cast<uint64_t>(c.qubits)));
        if (info(kind).nq == 2) {
            const auto q2 = static_cast<uint8_t>((q + 1 + r.bounded(static_cast<uint64_t>(c.qubits - 1))) % c.qubits);
            out.add(make_gate(kind, q, q2));
        } else {
            std::array<double,3> params{};
            for (int j = 0; j < info(kind).nparams; ++j) params[j] = angle(r, c.angle_mode);
            out.add(make_gate(kind, q, params));
        }
    }
    return out;
}
std::size_t planted(const GeneratorConfig& c) {
    return static_cast<std::size_t>(std::floor(c.prefix_fraction * static_cast<double>(c.gates)));
}
std::vector<bool> labels(const GeneratorConfig& c) {
    const auto count = static_cast<std::size_t>(std::floor(c.equivalent_fraction * static_cast<double>(c.pairs) + 0.5));
    std::vector<bool> out(c.pairs, false);
    for (std::size_t i = 0; i < count; ++i) out[i] = true;
    Random r{derive(c.seed, 4, 0)};
    for (std::size_t i = out.size(); i > 1; --i) {
        const auto j = static_cast<std::size_t>(r.bounded(i));
        const bool tmp = out[i-1]; out[i-1] = out[j]; out[j] = tmp;
    }
    return out;
}
std::string id(const char* prefix, std::size_t index) {
    std::ostringstream s; s << prefix << std::setw(6) << std::setfill('0') << index; return s.str();
}
bool equal_dense(const std::vector<cplx>& a, const std::vector<cplx>& b, bool ignore_phase, double tol) {
    cplx phase = 1;
    if (ignore_phase) {
        auto it = std::max_element(b.begin(), b.end(), [](cplx x, cplx y) { return std::norm(x) < std::norm(y); });
        const auto i = static_cast<std::size_t>(it - b.begin());
        phase = a[i] / b[i];
        if (std::abs(std::abs(phase) - 1) > tol) return false;
        phase /= std::abs(phase);
    }
    for (std::size_t i = 0; i < a.size(); ++i) if (std::abs(a[i] - phase * b[i]) > tol) return false;
    return true;
}
CircuitPair make_pair(const GeneratorConfig& c, std::size_t index, const Circuit& prefix, bool equivalent) {
    CircuitPair pair;
    auto& m = pair.metadata;
    m.case_id = id("case-", index);
    m.prefix_group_id = index / c.prefix_group_size;
    m.prefix_identity = id("prefix-", m.prefix_group_id);
    m.planted_prefix_gates = prefix.size();
    m.master_seed = c.seed;
    m.prefix_seed = derive(c.seed, 1, m.prefix_group_id);
    m.suffix_seed = derive(c.seed, 2, index);
    m.transform_seed = derive(c.seed, 3, index);
    pair.a = prefix;
    auto suffix = random_circuit(c, c.gates - prefix.size(), m.suffix_seed);
    pair.a.gates.insert(pair.a.gates.end(), suffix.gates.begin(), suffix.gates.end());
    pair.b = pair.a;
    m.transform = apply_equivalent_transform(pair.b, static_cast<TransformKind>(index % 7), prefix.size(), m.transform_seed);
    pair.expected = equivalent ? ExpectedRelation::Equivalent : ExpectedRelation::NonEquivalent;
    if (!equivalent) {
        Random r{derive(m.transform_seed, 5, 0)};
        m.negative_qubit = static_cast<int>(r.bounded(static_cast<uint64_t>(c.qubits)));
        pair.b.add(make_gate(GateKind::X, static_cast<uint8_t>(m.negative_qubit)));
    }
    m.sharing = measure_sharing(pair.a, pair.b);
    verify_pair(pair, c.verify);
    return pair;
}
} // namespace

const char* relation_name(ExpectedRelation r) { return r == ExpectedRelation::Equivalent ? "equivalent" : "non_equivalent"; }
const char* transform_name(TransformKind k) {
    switch (k) {
        case TransformKind::Identical: return "identical";
        case TransformKind::InsertInverse: return "insert_inverse";
        case TransformKind::XtoHZH: return "x_to_hzh";
        case TransformKind::StoTT: return "s_to_tt";
        case TransformKind::SwapToCX: return "swap_to_cx";
        case TransformKind::Commute: return "commute_disjoint";
        case TransformKind::GlobalPhase: return "global_phase";
    }
    throw std::invalid_argument("unknown transform kind");
}
void validate_config(const GeneratorConfig& c) {
    if (c.qubits < 1 || c.qubits > kMaxQubits) throw std::invalid_argument("qubits must be in [1,32]");
    // DAG uses uint32_t node indices. Leave room for inserted gates and the negative X.
    if (c.gates > UINT32_MAX - 4 || c.gates > std::vector<Gate>{}.max_size() - 4)
        throw std::invalid_argument("gate count exceeds IR/DAG capacity");
    if (!c.pairs || c.pairs > std::vector<CircuitPair>{}.max_size() || c.pairs > UINT32_MAX)
        throw std::invalid_argument("pairs must be positive and fit the batch capacity");
    if (!c.prefix_group_size) throw std::invalid_argument("prefix-group-size must be positive");
    for (double r : {c.prefix_fraction, c.equivalent_fraction})
        if (!std::isfinite(r) || r < 0 || r > 1) throw std::invalid_argument("fractions must be finite and in [0,1]");
    if (c.gate_set != "clifford" && c.gate_set != "clifford-t" && c.gate_set != "parameterized" && c.gate_set != "mixed")
        throw std::invalid_argument("unknown gate-set");
    if (c.angle_mode != "discrete" && c.angle_mode != "continuous" && c.angle_mode != "mixed")
        throw std::invalid_argument("unknown angle-mode");
    if (c.verify != "small" && c.verify != "none") throw std::invalid_argument("verify must be small or none");
}
void validate_circuit(const Circuit& c) {
    validate(c);
    if (c.size() > UINT32_MAX - 4) throw std::invalid_argument("circuit exceeds DAG capacity");
    for (const auto& g : c.gates) {
        const auto& gi = info(g.kind);
        for (int i = 0; i < 3; ++i) {
            if (!std::isfinite(g.params[i])) throw std::invalid_argument("non-finite gate parameter");
            if (i >= gi.nparams && g.params[i] != 0) throw std::invalid_argument("unused parameter must be zero");
        }
        if (g.kind == GateKind::U1Q || g.kind == GateKind::U2Q) {
            const int d = 1 << g.nq;
            for (int i = 0; i < d; ++i) for (int j = 0; j < d; ++j) {
                cplx v = 0;
                for (int k = 0; k < d; ++k) v += g.at(i,k) * std::conj(g.at(j,k));
                if (!std::isfinite(v.real()) || !std::isfinite(v.imag()) || std::abs(v - cplx(i == j ? 1 : 0)) > 1e-9)
                    throw std::invalid_argument("generic matrix is not unitary");
            }
        } else {
            const Gate expected = gi.nq == 1 ? make_gate(g.kind, g.q[0], g.params) : make_gate(g.kind, g.q[0], g.q[1]);
            if (!same_gate(g, expected, 1e-12)) throw std::invalid_argument("named gate matrix disagrees with parameters");
        }
    }
}
Circuit generate_circuit(const GeneratorConfig& c, uint64_t seed) {
    validate_config(c); return random_circuit(c, c.gates, seed);
}
Circuit inverse_circuit(const Circuit& c) {
    validate_circuit(c);
    Circuit out(c.n); out.global_phase = -c.global_phase;
    for (auto it = c.gates.rbegin(); it != c.gates.rend(); ++it) out.add(inverse(*it));
    return out;
}
TransformRecord apply_equivalent_transform(Circuit& c, TransformKind kind, std::size_t start, uint64_t seed) {
    validate_circuit(c);
    if (start > c.size()) throw std::invalid_argument("protected prefix exceeds circuit size");
    (void)transform_name(kind);
    Random r{seed};
    TransformRecord t; t.kind = kind; t.location = start;
    if (kind == TransformKind::Identical) { t.proof = "Copied circuit exactly."; return t; }
    if (kind == TransformKind::GlobalPhase) {
        t.phase_adjustment = pi / 7; c.global_phase += t.phase_adjustment;
        t.proof = "Only a circuit-wide scalar phase changes."; return t;
    }
    if (kind == TransformKind::XtoHZH || kind == TransformKind::StoTT || kind == TransformKind::SwapToCX) {
        const auto target = kind == TransformKind::XtoHZH ? GateKind::X : kind == TransformKind::StoTT ? GateKind::S : GateKind::SWAP;
        for (std::size_t i = start; i < c.size(); ++i) if (c.gates[i].kind == target) {
            const Gate g = c.gates[i]; t.location = i; t.removed = {g};
            if (target == GateKind::X) t.inserted = {make_gate(GateKind::H,g.q[0]),make_gate(GateKind::Z,g.q[0]),make_gate(GateKind::H,g.q[0])};
            if (target == GateKind::S) t.inserted = {make_gate(GateKind::T,g.q[0]),make_gate(GateKind::T,g.q[0])};
            if (target == GateKind::SWAP) t.inserted = {make_gate(GateKind::CX,g.q[0],g.q[1]),make_gate(GateKind::CX,g.q[1],g.q[0]),make_gate(GateKind::CX,g.q[0],g.q[1])};
            c.gates.erase(c.gates.begin() + static_cast<std::ptrdiff_t>(i));
            c.gates.insert(c.gates.begin() + static_cast<std::ptrdiff_t>(i), t.inserted.begin(), t.inserted.end());
            t.proof = target == GateKind::X ? "X = H Z H." : target == GateKind::S ? "S = T T." : "SWAP(a,b) = CX(a,b) CX(b,a) CX(a,b).";
            return t;
        }
    }
    if (kind == TransformKind::Commute) {
        for (std::size_t i = start; i + 1 < c.size(); ++i) if (!overlaps(c.gates[i],c.gates[i+1])) {
            t.location = i; t.removed = {c.gates[i],c.gates[i+1]};
            std::swap(c.gates[i],c.gates[i+1]); t.inserted = {c.gates[i],c.gates[i+1]};
            t.proof = "Adjacent gates on disjoint qubits commute."; return t;
        }
    }
    // Empty suffixes retain the entire base prefix. Ineligible rewrites fall back to G G^-1.
    const bool fallback = kind != TransformKind::InsertInverse;
    t.kind = TransformKind::InsertInverse;
    t.location = start + static_cast<std::size_t>(r.bounded(c.size() - start + 1));
    const auto q = static_cast<uint8_t>(r.bounded(static_cast<uint64_t>(c.n)));
    const Gate g = c.gates.empty() ? make_gate(GateKind::H,q) : c.gates[static_cast<std::size_t>(r.bounded(c.size()))];
    t.inserted = {g,inverse(g)};
    c.gates.insert(c.gates.begin() + static_cast<std::ptrdiff_t>(t.location),t.inserted.begin(),t.inserted.end());
    t.proof = fallback ? "No eligible suffix rewrite; inserted G G^-1 = I." : "Inserted G G^-1 = I.";
    return t;
}
PrefixMetrics common_prefix(const Circuit& a, const Circuit& b, double tol) {
    if (!std::isfinite(tol) || tol < 0) throw std::invalid_argument("invalid prefix tolerance");
    PrefixMetrics m; m.a_gates = a.size(); m.b_gates = b.size();
    const auto n = std::min(a.size(),b.size());
    if (a.n == b.n) while (m.gates < n && same_gate(a.gates[m.gates],b.gates[m.gates],tol)) ++m.gates;
    m.fraction = n ? static_cast<double>(m.gates)/static_cast<double>(n) : (a.size() == b.size() ? 1.0 : 0.0);
    return m;
}
Circuit normalized_copy(const Circuit& c) {
    validate_circuit(c); CircuitDAG dag(c); normalize(dag);
    Circuit out = dag.to_circuit(); validate_circuit(out); return out;
}
SharingMetrics measure_sharing(const Circuit& a, const Circuit& b) {
    validate_circuit(a); validate_circuit(b);
    return {common_prefix(a,b),common_prefix(normalized_copy(a),normalized_copy(b))};
}
std::vector<cplx> reference_unitary(const Circuit& c) {
    validate_circuit(c);
    if (c.n > kDenseQubitCap) throw std::invalid_argument("dense verification limited to 5 qubits");
    const std::size_t d = std::size_t{1} << c.n;
    std::vector<cplx> u(d*d);
    for (std::size_t col = 0; col < d; ++col) {
        cplx* v = u.data() + col*d; v[col] = std::polar(1.0,c.global_phase);
        for (const auto& g : c.gates) {
            const std::size_t b0 = std::size_t{1} << g.q[0];
            const std::size_t b1 = g.nq == 2 ? std::size_t{1} << g.q[1] : 0;
            const int local = 1 << g.nq;
            for (std::size_t base = 0; base < d; ++base) {
                if (base & (b0 | b1)) continue;
                const std::size_t ix[4] = {base,base|b0,base|b1,base|b0|b1};
                cplx old[4]{}, out[4]{};
                for (int j = 0; j < local; ++j) old[j] = v[ix[j]];
                for (int i = 0; i < local; ++i) for (int j = 0; j < local; ++j) out[i] += g.at(i,j)*old[j];
                for (int i = 0; i < local; ++i) v[ix[i]] = out[i];
            }
        }
    }
    return u;
}
bool equivalent_up_to_phase(const Circuit& a, const Circuit& b, double tol) {
    if (a.n != b.n) return false;
    if (!std::isfinite(tol) || tol <= 0) throw std::invalid_argument("invalid verification tolerance");
    return equal_dense(reference_unitary(a),reference_unitary(b),true,tol);
}
void verify_pair(CircuitPair& p, const std::string& mode) {
    validate_circuit(p.a); validate_circuit(p.b);
    if (p.a.n != p.b.n) throw std::invalid_argument("pair qubit counts differ");
    if (mode == "none") { p.metadata.verification = "skipped_by_request"; return; }
    if (mode != "small") throw std::invalid_argument("unknown verification mode");
    if (p.a.n > kDenseQubitCap) { p.metadata.verification = "skipped_size_limit"; return; }
    const auto a = reference_unitary(p.a), b = reference_unitary(p.b);
    if (equal_dense(a,b,true,kVerificationTolerance) != (p.expected == ExpectedRelation::Equivalent))
        throw std::runtime_error("label verification failed for " + p.metadata.case_id);
    if (!equal_dense(a,reference_unitary(normalized_copy(p.a)),false,kVerificationTolerance) ||
        !equal_dense(b,reference_unitary(normalized_copy(p.b)),false,kVerificationTolerance))
        throw std::runtime_error("normalization changed unitary for " + p.metadata.case_id);
    p.metadata.verification = "passed";
}
CircuitPair generate_pair(const GeneratorConfig& c, std::size_t index) {
    validate_config(c);
    if (index >= c.pairs) throw std::invalid_argument("case index outside batch");
    const auto prefix = random_circuit(c,planted(c),derive(c.seed,1,index/c.prefix_group_size));
    return make_pair(c,index,prefix,labels(c)[index]);
}
std::vector<CircuitPair> generate_batch(const GeneratorConfig& c) {
    validate_config(c);
    const auto expected = labels(c);
    std::vector<CircuitPair> out; out.reserve(c.pairs);
    Circuit prefix(c.qubits);
    for (std::size_t i = 0; i < c.pairs; ++i) {
        if (i % c.prefix_group_size == 0) prefix = random_circuit(c,planted(c),derive(c.seed,1,i/c.prefix_group_size));
        out.push_back(make_pair(c,i,prefix,expected[i]));
    }
    return out;
}
} // namespace qcanopy::testgen
