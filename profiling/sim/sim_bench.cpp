// sim_bench: one profiling check of qcanopy::Simulator.
//
// Pair mode (one labeled pair, run by sim_profile.py in a fresh process):
//   sim_bench --a A.json --b B.json [--strategy flat|dag|normalized]
//             [--repeats 3] [--probes 3] [--probe-seed 759] [--tol 1e-8]
//             [--reference-max-qubits 12] [--reference-max-work 400000000]
//             [--no-warmup]
//
// Kernel-scan mode (per-gate cost by gate kind and target qubit):
//   sim_bench --kernel-scan --qubits N [--repeats 3]
//
// Output is one JSON object on stdout. Exit 0 on success, 2 on bad arguments,
// 3 on a runtime error (the JSON then carries "error").

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#elif defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "qcanopy/sim.hpp"
#include "sim_bench_adapter.hpp"

using qcanopy::Circuit;
using qcanopy::Simulator;
using qcanopy::StateVector;
using Clock = std::chrono::steady_clock;
using cd = std::complex<double>;

namespace {

// ---------------------------------------------------------------- utilities

double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

double median(std::vector<double> v) {
    if (v.empty()) return std::nan("");
    std::sort(v.begin(), v.end());
    size_t m = v.size() / 2;
    return v.size() % 2 ? v[m] : 0.5 * (v[m - 1] + v[m]);
}

// SplitMix64, same constants as the generator (docs/generator.md).
struct SplitMix64 {
    uint64_t state;
    static uint64_t mix(uint64_t z) {
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    uint64_t next() { state += 0x9e3779b97f4a7c15ULL; return mix(state); }
    double uniform() { return double(next() >> 11) * 0x1.0p-53; }
};

// ------------------------------------------------------------- memory (RSS)

struct Rss { double current_mb = std::nan(""), peak_mb = std::nan(""); const char* source = "none"; };

#if defined(__linux__)
double read_status_kb(const char* key) {
    std::ifstream f("/proc/self/status");
    std::string line;
    size_t klen = std::strlen(key);
    while (std::getline(f, line)) {
        if (line.compare(0, klen, key) == 0) return std::atof(line.c_str() + klen + 1);
    }
    return std::nan("");
}
#endif

void reset_peak_rss() {
#if defined(__linux__)
    // Writing 5 resets VmHWM to the current RSS (Linux >= 4.0).
    std::ofstream f("/proc/self/clear_refs");
    if (f) f << "5";
#endif
}

Rss read_rss() {
    Rss r;
#if defined(__linux__)
    r.current_mb = read_status_kb("VmRSS:") / 1024.0;
    r.peak_mb = read_status_kb("VmHWM:") / 1024.0;
    r.source = "VmHWM";
#elif defined(_WIN32)
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        r.current_mb = double(pmc.WorkingSetSize) / (1024.0 * 1024.0);
        r.peak_mb = double(pmc.PeakWorkingSetSize) / (1024.0 * 1024.0);
        r.source = "PeakWorkingSetSize";
    }
#elif defined(__APPLE__)
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) == 0) {
        r.peak_mb = double(ru.ru_maxrss) / (1024.0 * 1024.0);  // bytes on macOS
        r.source = "ru_maxrss";
    }
#endif
    return r;
}

// --------------------------------------------------------------- JSON output

std::string jnum(double v) {
    if (!std::isfinite(v)) return "null";
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.10g", v);
    return buf;
}

std::string jstr(const std::string& s) {
    std::string o = "\"";
    for (char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\r': o += "\\r"; break;
            case '\t': o += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char b[8];
                    std::snprintf(b, sizeof b, "\\u%04x", c);
                    o += b;
                } else {
                    o += c;
                }
        }
    }
    return o + "\"";
}

std::string jarr(const std::vector<double>& v) {
    std::string o = "[";
    for (size_t i = 0; i < v.size(); ++i) o += (i ? "," : "") + jnum(v[i]);
    return o + "]";
}

// ------------------------------------------------------------ circuit tools

Circuit load_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return qcanopy::simbench::load_ir(f);
}

Circuit load_string(const std::string& json) {
    std::istringstream in(json);
    return qcanopy::simbench::load_ir(in);
}

// A layer of random U gates, one per qubit, as IR JSON. Prepending the same
// layer to A and B keeps equivalence (A P ~ B P iff A ~ B) and turns the
// fixed |0...0> input into a random product state.
Circuit probe_layer(int n, uint64_t seed) {
    SplitMix64 rng{seed};
    const double pi = 3.14159265358979323846;
    std::ostringstream js;
    js.precision(17);
    js << "{\"schema_version\":1,\"n\":" << n << ",\"global_phase\":0,\"gates\":[";
    for (int q = 0; q < n; ++q) {
        double theta = std::acos(1.0 - 2.0 * rng.uniform());  // uniform on the sphere
        double phi = 2.0 * pi * rng.uniform();
        double lam = 2.0 * pi * rng.uniform();
        js << (q ? "," : "") << "{\"name\":\"u\",\"qubits\":[" << q << "],\"params\":["
           << theta << "," << phi << "," << lam << "]}";
    }
    js << "]}";
    return load_string(js.str());
}

Circuit with_prefix(const Circuit& prefix, const Circuit& c) {
    Circuit out = prefix;
    out.gates.insert(out.gates.end(), c.gates.begin(), c.gates.end());
    out.global_phase = prefix.global_phase + c.global_phase;
    return out;
}

// Independent reference: a deliberately naive dense application that loops
// over every basis index and does not reuse sim.cpp's index-insertion trick.
// Local index convention: (bit(q[1]) << 1) | bit(q[0]); matrices row-major.
std::vector<cd> reference_state(const Circuit& c) {
    const size_t dim = size_t(1) << c.n;
    std::vector<cd> a(dim, cd(0, 0)), out(dim);
    a[0] = cd(1, 0);
    for (const auto& g : c.gates) {
        std::fill(out.begin(), out.end(), cd(0, 0));
        if (g.nq == 1) {
            const size_t b0 = size_t(1) << g.q[0];
            for (size_t i = 0; i < dim; ++i) {
                if (a[i] == cd(0, 0)) continue;
                const size_t col = (i & b0) ? 1 : 0;
                const size_t base = i & ~b0;
                out[base] += g.m[0 * 2 + col] * a[i];
                out[base | b0] += g.m[1 * 2 + col] * a[i];
            }
        } else if (g.nq == 2) {
            const size_t b0 = size_t(1) << g.q[0], b1 = size_t(1) << g.q[1];
            for (size_t i = 0; i < dim; ++i) {
                if (a[i] == cd(0, 0)) continue;
                const size_t col = ((i & b1) ? 2 : 0) | ((i & b0) ? 1 : 0);
                const size_t base = i & ~(b0 | b1);
                for (size_t row = 0; row < 4; ++row) {
                    const size_t j = base | ((row & 1) ? b0 : 0) | ((row & 2) ? b1 : 0);
                    out[j] += g.m[row * 4 + col] * a[i];
                }
            }
        } else {
            throw std::runtime_error("reference supports 1- and 2-qubit gates only");
        }
        a.swap(out);
    }
    const cd ph = std::polar(1.0, c.global_phase);
    for (auto& x : a) x *= ph;
    return a;
}

// ------------------------------------------------------------- strategies

enum class Strategy { Flat, Dag, Normalized };

// A circuit prepared once (untimed in sim_s, timed as prep_s), then executed.
struct Prepared {
    Strategy strategy;
    Circuit circuit;  // flat / normalized
    std::vector<qcanopy::CircuitDAG> dag;  // 0 or 1 element (no default ctor assumed)
    double prep_s = 0.0;
    size_t gates = 0;  // gates the simulator actually applies
};

Prepared prepare(Strategy s, const Circuit& c) {
    Prepared p{s, c, {}, 0.0, c.gates.size()};
    auto t0 = Clock::now();
    if (s == Strategy::Normalized) {
        p.circuit = qcanopy::simbench::normalized(c);
        p.prep_s = seconds_since(t0);
        p.gates = p.circuit.gates.size();
    } else if (s == Strategy::Dag) {
        p.dag.push_back(qcanopy::simbench::make_dag(c));
        p.prep_s = seconds_since(t0);
    }
    return p;
}

StateVector run(const Simulator& sim, const Prepared& p) {
    if (p.strategy == Strategy::Dag) return sim.execute(p.dag.front());
    return sim.execute(p.circuit);
}

struct Overlap { cd inner; double norm_a, norm_b; };

Overlap overlap(const StateVector& a, const StateVector& b) {
    cd s(0, 0);
    double na = 0, nb = 0;
    for (size_t i = 0; i < a.amplitudes.size(); ++i) {
        s += std::conj(a.amplitudes[i]) * b.amplitudes[i];
        na += std::norm(a.amplitudes[i]);
        nb += std::norm(b.amplitudes[i]);
    }
    return {s, std::sqrt(na), std::sqrt(nb)};
}

// ------------------------------------------------------------- arguments

struct Args {
    std::string a, b, strategy = "flat";
    int repeats = 3, probes = 3, reference_max_qubits = 12, qubits = -1;
    uint64_t probe_seed = 759;
    double tol = 1e-8, reference_max_work = 4e8;
    bool warmup = true, kernel_scan = false;
};

[[noreturn]] void usage(const std::string& msg) {
    std::cout << "{\"error\":" << jstr("usage: " + msg) << "}\n";
    std::exit(2);
}

Args parse(int argc, char** argv) {
    Args a;
    auto need = [&](int& i) -> std::string {
        if (i + 1 >= argc) usage(std::string("missing value for ") + argv[i]);
        return argv[++i];
    };
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        try {
            if (k == "--a") a.a = need(i);
            else if (k == "--b") a.b = need(i);
            else if (k == "--strategy") a.strategy = need(i);
            else if (k == "--repeats") a.repeats = std::stoi(need(i));
            else if (k == "--probes") a.probes = std::stoi(need(i));
            else if (k == "--probe-seed") a.probe_seed = std::stoull(need(i));
            else if (k == "--tol") a.tol = std::stod(need(i));
            else if (k == "--reference-max-qubits") a.reference_max_qubits = std::stoi(need(i));
            else if (k == "--reference-max-work") a.reference_max_work = std::stod(need(i));
            else if (k == "--no-warmup") a.warmup = false;
            else if (k == "--kernel-scan") a.kernel_scan = true;
            else if (k == "--qubits") a.qubits = std::stoi(need(i));
            else usage("unknown option " + k);
        } catch (const std::logic_error&) {
            usage("bad value for " + k);
        }
    }
    if (a.repeats < 1) usage("--repeats must be >= 1");
    if (a.kernel_scan) {
        if (a.qubits < 1 || a.qubits > 34) usage("--kernel-scan needs --qubits 1..34");
    } else {
        if (a.a.empty() || a.b.empty()) usage("--a and --b are required");
        if (a.probes < 1) usage("--probes must be >= 1");
        if (a.strategy != "flat" && a.strategy != "dag" && a.strategy != "normalized")
            usage("--strategy must be flat, dag or normalized");
    }
    return a;
}

void warm_up(const Simulator& sim) {
    Circuit w = load_string(
        "{\"schema_version\":1,\"n\":2,\"global_phase\":0.1,\"gates\":["
        "{\"name\":\"h\",\"qubits\":[0],\"params\":[]},"
        "{\"name\":\"cx\",\"qubits\":[0,1],\"params\":[]}]}");
    volatile double sink = sim.execute(w).probability(3);
    (void)sink;
}

// ------------------------------------------------------------- pair mode

int pair_mode(const Args& args) {
    const Simulator sim;
    const Strategy strategy = args.strategy == "dag" ? Strategy::Dag
                            : args.strategy == "normalized" ? Strategy::Normalized
                            : Strategy::Flat;
    if (args.warmup) warm_up(sim);

    const Rss before = read_rss();
    reset_peak_rss();

    auto t0 = Clock::now();
    Circuit A = load_file(args.a);
    Circuit B = load_file(args.b);
    const double load_s = seconds_since(t0);
    if (A.n != B.n) throw std::runtime_error("A and B have different qubit counts");
    const int n = A.n;

    // Probe 0: the |0...0> input exactly as the simulator sees it. Timed.
    Prepared pa = prepare(strategy, A), pb = prepare(strategy, B);
    std::vector<double> sim_a, sim_b;
    StateVector sa(1), sb(1);
    for (int r = 0; r < args.repeats; ++r) {
        auto t = Clock::now();
        sa = run(sim, pa);
        sim_a.push_back(seconds_since(t));
        t = Clock::now();
        sb = run(sim, pb);
        sim_b.push_back(seconds_since(t));
    }

    // Reference check of probe 0 (exact amplitudes, including global phase).
    double ref_err_a = std::nan(""), ref_err_b = std::nan("");
    bool ref_ran = false;
    const double ref_work = double(A.gates.size() + B.gates.size()) * std::ldexp(1.0, n);
    if (n <= args.reference_max_qubits && ref_work <= args.reference_max_work) {
        ref_ran = true;
        auto cmp = [](const std::vector<cd>& ref, const StateVector& s) {
            double e = 0;
            for (size_t i = 0; i < ref.size(); ++i) e = std::max(e, std::abs(ref[i] - s.amplitudes[i]));
            return e;
        };
        ref_err_a = cmp(reference_state(A), sa);
        ref_err_b = cmp(reference_state(B), sb);
    }

    // Probes: overlap <A psi | B psi>. Equivalent up to one global phase means
    // |overlap| = 1 for every probe AND the same phase on every probe.
    struct ProbeRow { double fidelity, phase, norm_a, norm_b; };
    std::vector<ProbeRow> rows;
    auto record = [&](const StateVector& x, const StateVector& y) {
        Overlap o = overlap(x, y);
        rows.push_back({std::abs(o.inner), std::arg(o.inner), o.norm_a, o.norm_b});
    };
    record(sa, sb);
    {   // free probe-0 states before allocating more
        StateVector e(1);
        sa = e; sb = e;
    }
    for (int k = 1; k < args.probes; ++k) {
        Circuit P = probe_layer(n, SplitMix64::mix(args.probe_seed + uint64_t(k)));
        Prepared qa = prepare(strategy, with_prefix(P, A));
        Prepared qb = prepare(strategy, with_prefix(P, B));
        StateVector xa = run(sim, qa);
        StateVector xb = run(sim, qb);
        record(xa, xb);
    }

    const double pi = 3.14159265358979323846;
    double min_fid = 2, max_norm_err = 0, phase_spread = 0;
    for (const auto& r : rows) {
        min_fid = std::min(min_fid, r.fidelity);
        max_norm_err = std::max({max_norm_err, std::abs(r.norm_a - 1), std::abs(r.norm_b - 1)});
        double d = std::remainder(r.phase - rows.front().phase, 2 * pi);
        phase_spread = std::max(phase_spread, std::abs(d));
    }
    const double phase_tol = std::max(1e-6, std::sqrt(args.tol));
    const bool equivalent = min_fid >= 1 - args.tol && phase_spread <= phase_tol;

    const Rss after = read_rss();

    std::ostringstream o;
    o << "{\"mode\":\"pair\",\"n\":" << n
      << ",\"strategy\":" << jstr(args.strategy)
      << ",\"gates_a\":" << A.gates.size() << ",\"gates_b\":" << B.gates.size()
      << ",\"simulated_gates_a\":" << pa.gates << ",\"simulated_gates_b\":" << pb.gates
      << ",\"state_bytes\":" << jnum(16.0 * std::ldexp(1.0, n))
      << ",\"load_s\":" << jnum(load_s)
      << ",\"prep_a_s\":" << jnum(pa.prep_s) << ",\"prep_b_s\":" << jnum(pb.prep_s)
      << ",\"sim_a_s\":" << jarr(sim_a) << ",\"sim_b_s\":" << jarr(sim_b)
      << ",\"sim_a_median_s\":" << jnum(median(sim_a)) << ",\"sim_b_median_s\":" << jnum(median(sim_b))
      << ",\"probes\":[";
    for (size_t i = 0; i < rows.size(); ++i) {
        o << (i ? "," : "") << "{\"fidelity\":" << jnum(rows[i].fidelity)
          << ",\"phase\":" << jnum(rows[i].phase) << ",\"norm_a\":" << jnum(rows[i].norm_a)
          << ",\"norm_b\":" << jnum(rows[i].norm_b) << "}";
    }
    o << "],\"min_fidelity\":" << jnum(min_fid) << ",\"phase_spread\":" << jnum(phase_spread)
      << ",\"max_norm_err\":" << jnum(max_norm_err)
      << ",\"verdict\":" << jstr(equivalent ? "equivalent" : "not_equivalent")
      << ",\"reference\":" << (ref_ran ? "true" : "false")
      << ",\"ref_err_a\":" << jnum(ref_err_a) << ",\"ref_err_b\":" << jnum(ref_err_b)
      << ",\"rss_before_mb\":" << jnum(before.current_mb)
      << ",\"peak_rss_mb\":" << jnum(after.peak_mb)
      << ",\"rss_source\":" << jstr(after.source) << "}";
    std::cout << o.str() << "\n";
    return 0;
}

// ------------------------------------------------------------- kernel scan

// Times circuits made of G copies of one gate. Per-gate cost is
// (t(G gates) - t(empty circuit)) / G, which removes state allocation.
int kernel_scan_mode(const Args& args) {
    const Simulator sim;
    const int n = args.qubits;
    const double dim = std::ldexp(1.0, n);
    // Keep each timed circuit around 2^23 amplitude updates (bounded 2..4096 gates).
    const int G = int(std::clamp(std::ldexp(1.0, 23) / dim, 2.0, 4096.0));
    if (args.warmup) warm_up(sim);

    auto header = [&](std::ostringstream& js) {
        js << "{\"schema_version\":1,\"n\":" << n << ",\"global_phase\":0,\"gates\":[";
    };
    auto time_circuit = [&](const Circuit& c) {
        std::vector<double> t;
        for (int r = 0; r < args.repeats; ++r) {
            auto t0 = Clock::now();
            StateVector s = sim.execute(c);
            t.push_back(seconds_since(t0));
            volatile double sink = s.probability(0);
            (void)sink;
        }
        return median(t);
    };

    std::ostringstream ej;
    header(ej);
    ej << "]}";
    const double empty_s = time_circuit(load_string(ej.str()));

    struct Spec { const char* kind; int q0, q1; };
    std::vector<Spec> specs;
    // Dense kernels (h, cx) at every target; diagonal / permutation kernels
    // (z, cz, swap) at the lowest, middle and highest target for comparison.
    auto few = [n](int q) { return q == 0 || q == 1 || q == n / 2 || q == n - 1; };
    for (int q = 0; q < n; ++q) {
        specs.push_back({"h", q, -1});
        if (few(q)) specs.push_back({"z", q, -1});
    }
    for (int q = 1; q < n; ++q) {
        specs.push_back({"cx", 0, q});
        if (few(q)) { specs.push_back({"cz", 0, q}); specs.push_back({"swap", 0, q}); }
    }

    std::ostringstream o;
    o << "{\"mode\":\"kernel_scan\",\"n\":" << n << ",\"gates_per_circuit\":" << G
      << ",\"empty_circuit_s\":" << jnum(empty_s) << ",\"rows\":[";
    for (size_t i = 0; i < specs.size(); ++i) {
        const Spec& s = specs[i];
        std::ostringstream js;
        header(js);
        for (int g = 0; g < G; ++g) {
            js << (g ? "," : "") << "{\"name\":\"" << s.kind << "\",\"qubits\":[" << s.q0;
            if (s.q1 >= 0) js << "," << s.q1;
            js << "],\"params\":[]}";
        }
        js << "]}";
        const double t = time_circuit(load_string(js.str()));
        const double per_gate = std::max(0.0, t - empty_s) / G;
        o << (i ? "," : "") << "{\"kind\":" << jstr(s.kind) << ",\"nq\":" << (s.q1 >= 0 ? 2 : 1)
          << ",\"target\":" << (s.q1 >= 0 ? s.q1 : s.q0)
          << ",\"ns_per_gate\":" << jnum(per_gate * 1e9)
          << ",\"ns_per_amp\":" << jnum(per_gate * 1e9 / dim)
          << ",\"gb_s\":" << jnum(per_gate > 0 ? 32.0 * dim / per_gate / 1e9 : std::nan("")) << "}";
    }
    o << "]}";
    std::cout << o.str() << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    Args args = parse(argc, argv);
    try {
        return args.kernel_scan ? kernel_scan_mode(args) : pair_mode(args);
    } catch (const std::bad_alloc&) {
        std::cout << "{\"error\":\"out of memory (std::bad_alloc)\"}\n";
    } catch (const std::exception& e) {
        std::cout << "{\"error\":" << jstr(e.what()) << "}\n";
    }
    return 3;
}
