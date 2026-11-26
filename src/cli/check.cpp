#include "check.hpp"

#include "qcanopy/dag.hpp"
#include "qcanopy/normalize.hpp"
#include "qcanopy/parser.hpp"
#include "qcanopy/peel.hpp"
#include "qcanopy/stimulus.hpp"

#include <charconv>
#include <cmath>
#include <complex>
#include <exception>
#include <iomanip>
#include <iostream>
#include <locale>
#include <set>
#include <sstream>

namespace qcanopy::cli {

    const char* const check_help =
        "qcanopy check A.qasm B.qasm [options]\n"
        "Check whether two OpenQASM 3.0 circuits are equal up to a global phase.\n"
        "\n"
        "  --stimulus LIST   Comma-separated input states to simulate (default haar,product):\n"
        "                      haar       Haar-random states (--count of them)\n"
        "                      product    random product states (--count of them)\n"
        "                      basis      random computational basis states (--count of them)\n"
        "                      all-basis  every basis state; passing them is a proof (at most 12 qubits)\n"
        "  --count N         States per random kind (default 8).\n"
        "  --seed UINT64     Seed for the random kinds (default 759).\n"
        "  --tol X           Tolerance on the residual norm and phase (default 1e-9).\n"
        "  --no-rewrite      Skip the normalize + peel pre-pass and simulate the full circuits.\n"
        "                    The pre-pass matches gates to about 1e-9, so use this with a smaller --tol.\n"
        "  --json            Print the result as one JSON object.\n"
        "  --quiet           Print nothing; report the verdict by exit status only.\n"
        "  --help            Show this help.\n"
        "\n"
        "Exit status: 0 equivalent, 1 not equivalent, 2 usage, parse or I/O error.\n"
        "A non-equivalent verdict is definitive. An equivalent verdict from random stimuli means\n"
        "no stimulus told the circuits apart; from all-basis or the pre-pass it is a proof.\n";

    namespace {
        const std::set<std::string> kStimulusKinds{"haar", "product", "basis", "all-basis"};

        uint64_t parse_uint(const std::string& key, const std::string& s) {
            uint64_t out = 0;
            const auto r = std::from_chars(s.data(), s.data() + s.size(), out);
            if (s.empty() || r.ec != std::errc{} || r.ptr != s.data() + s.size())
                throw UsageError(key + " expects a non-negative integer, got '" + s + "'");
            return out;
        }

        double parse_real(const std::string& key, const std::string& s) {
            std::istringstream in(s);
            in.imbue(std::locale::classic());
            double out = 0.0;
            in >> std::noskipws >> out;
            if (!in || !in.eof() || !std::isfinite(out))
                throw UsageError(key + " expects a finite number, got '" + s + "'");
            return out;
        }

        std::vector<std::string> parse_stimuli(const std::string& s) {
            std::vector<std::string> out;
            std::set<std::string> seen;
            std::size_t start = 0;
            while (true) {
                const std::size_t comma = s.find(',', start);
                const std::string item = s.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
                if (!kStimulusKinds.count(item))
                    throw UsageError("unknown stimulus '" + item + "' (expected haar, product, basis or all-basis)");
                if (seen.insert(item).second) out.push_back(item);
                if (comma == std::string::npos) break;
                start = comma + 1;
            }
            return out;
        }

        const char* kind_name(StimulusKind k) {
            switch (k) {
                case StimulusKind::Haar: return "haar";
                case StimulusKind::Product: return "product";
                case StimulusKind::Basis: return "basis";
                case StimulusKind::Custom: break;
            }
            return "custom";
        }

        const char* mismatch_name(Mismatch m) {
            switch (m) {
                case Mismatch::State: return "state";
                case Mismatch::Phase: return "phase";
                case Mismatch::None: break;
            }
            return "none";
        }

        std::string json_string(const std::string& s) {
            std::ostringstream o;
            o << '"';
            for (const unsigned char c : s) {
                if (c == '"' || c == '\\') o << '\\' << c;
                else if (c < 0x20) o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
                else o << c;
            }
            o << '"';
            return o.str();
        }

        bool has_all_basis(const CheckOptions& opts) {
            for (const auto& s : opts.stimuli) if (s == "all-basis") return true;
            return false;
        }

        void print_text(std::ostream& out, const CheckOptions& opts, const CheckOutcome& r) {
            const EquivResult& e = r.equiv;
            out << std::setprecision(6);
            if (e.equivalent) {
                double phase = std::arg(e.global_phase);
                if (std::abs(phase) < 1e-12) phase = 0.0;  // hide rounding noise like -2.9e-17
                out << "equivalent up to global phase e^{i*" << phase << "}";
                if (r.proved_by_rewrite) {
                    out << " (proved by rewriting: no gates left after normalize + peel)\n";
                    return;
                }
                out << " (" << e.stimuli_checked << " stimuli, max deviation " << std::setprecision(2)
                    << e.deviation << "; " << (has_all_basis(opts) ? "all basis states: proof" : "not a proof")
                    << ")\n";
            } else {
                out << "NOT equivalent: stimulus #" << e.failing_index << " (" << kind_name(e.failing_kind) << "): "
                    << (e.mismatch == Mismatch::State ? "output is not a phase multiple of the input"
                                                      : "output phase differs from stimulus #0's")
                    << " (deviation " << std::setprecision(2) << e.deviation << ")\n";
            }
            if (opts.rewrite)
                out << "  gates: " << r.gates_a << " vs " << r.gates_b << ", after rewriting " << r.residual_a
                    << " vs " << r.residual_b << "\n";
        }

        void print_json(std::ostream& out, const CheckOptions& opts, const CheckOutcome& r) {
            const EquivResult& e = r.equiv;
            out << std::setprecision(17);
            out << "{\"a\":" << json_string(opts.file_a) << ",\"b\":" << json_string(opts.file_b)
                << ",\"qubits\":" << r.qubits << ",\"equivalent\":" << (e.equivalent ? "true" : "false")
                << ",\"method\":\"" << (r.proved_by_rewrite ? "rewrite" : "simulation") << "\""
                << ",\"gates\":[" << r.gates_a << "," << r.gates_b << "]"
                << ",\"residual_gates\":[" << r.residual_a << "," << r.residual_b << "]"
                << ",\"stimuli\":" << r.stimuli << ",\"stimuli_checked\":" << e.stimuli_checked
                << ",\"deviation\":" << e.deviation;
            if (e.equivalent) out << ",\"global_phase\":" << std::arg(e.global_phase);
            else
                out << ",\"witness\":{\"index\":" << e.failing_index << ",\"kind\":\"" << kind_name(e.failing_kind)
                    << "\",\"mismatch\":\"" << mismatch_name(e.mismatch) << "\"}";
            out << "}\n";
        }
    } // namespace

    CheckOptions parse_check_args(int argc, const char* const* argv) {
        CheckOptions opts;
        std::vector<std::string> files;
        std::set<std::string> seen;
        for (int i = 1; i < argc; ++i) {
            const std::string key = argv[i];
            if (key == "--help") { opts.help = true; return opts; }
            if (key.rfind("--", 0) != 0) { files.push_back(key); continue; }
            if (!seen.insert(key).second) throw UsageError("duplicate option " + key);

            if (key == "--no-rewrite") { opts.rewrite = false; continue; }
            if (key == "--json") { opts.json = true; continue; }
            if (key == "--quiet") { opts.quiet = true; continue; }

            if (key != "--stimulus" && key != "--count" && key != "--seed" && key != "--tol")
                throw UsageError("unknown option " + key);
            if (i + 1 == argc) throw UsageError("missing value for " + key);
            const std::string v = argv[++i];
            if (key == "--stimulus") opts.stimuli = parse_stimuli(v);
            else if (key == "--count") opts.count = static_cast<std::size_t>(parse_uint(key, v));
            else if (key == "--seed") opts.seed = parse_uint(key, v);
            else opts.tol = parse_real(key, v);
        }
        if (files.size() != 2)
            throw UsageError("expected two circuit files, got " + std::to_string(files.size()));
        if (opts.json && opts.quiet) throw UsageError("--json and --quiet cannot be combined");
        if (opts.count == 0) throw UsageError("--count must be at least 1");
        if (!(opts.tol > 0.0)) throw UsageError("--tol must be positive");
        opts.file_a = files[0];
        opts.file_b = files[1];
        return opts;
    }

    Stimulus build_stimulus(int n, const CheckOptions& opts) {
        Stimulus stim(n);
        for (const auto& kind : opts.stimuli) {
            if (kind == "haar") stim.append(haar_random_states(n, opts.count, opts.seed));
            else if (kind == "product") stim.append(random_product_states(n, opts.count, opts.seed));
            else if (kind == "basis") stim.append(random_basis_states(n, opts.count, opts.seed));
            else if (kind == "all-basis") {
                if (n > kMaxAllBasisQubits)
                    throw UsageError("all-basis is limited to " + std::to_string(kMaxAllBasisQubits) +
                                     " qubits; these circuits have " + std::to_string(n));
                stim.append(all_basis_states(n));
            } else {
                throw UsageError("unknown stimulus '" + kind + "'");
            }
        }
        return stim;
    }

    CheckOutcome check_circuits(const Circuit& a, const Circuit& b, const CheckOptions& opts) {
        if (a.n != b.n) throw std::invalid_argument("circuits have different qubit counts");

        CheckOutcome r;
        r.qubits = a.n;
        r.gates_a = a.size();
        r.gates_b = b.size();

        // The miter is B^-1 A, so B is peel's reference R and A its candidate C: peeled
        // phases land in A's global_phase and the residual miter keeps phase(A) - phase(B).
        Circuit a_mid = a, b_mid = b;
        if (opts.rewrite) {
            CircuitDAG da(a), db(b);
            normalize(da);
            normalize(db);
            // peel() already removes the maximal prefix and suffix, so this runs at most twice
            // today. It is the hook for a future rewriter, whose changes can expose new matches.
            while (peel(/*R=*/db, /*C=*/da) > 0) {}
            a_mid = da.to_circuit();
            b_mid = db.to_circuit();
        }
        r.residual_a = a_mid.size();
        r.residual_b = b_mid.size();

        if (a_mid.gates.empty() && b_mid.gates.empty()) {
            r.proved_by_rewrite = opts.rewrite;
            r.equiv.equivalent = true;
            r.equiv.global_phase = std::polar(1.0, a_mid.global_phase - b_mid.global_phase);
            return r;
        }

        const Stimulus stim = build_stimulus(a.n, opts);
        r.stimuli = stim.size();
        r.equiv = EquivalenceChecker(opts.tol).check(a_mid, b_mid, stim);
        return r;
    }

    int run_check(const CheckOptions& opts, std::ostream& out, std::ostream& err) {
        try {
            const Circuit a = parse_qasm_file(opts.file_a);
            const Circuit b = parse_qasm_file(opts.file_b);
            if (a.n != b.n) {
                err << "qcanopy check: " << opts.file_a << " has " << a.n << " qubits but " << opts.file_b
                    << " has " << b.n << "\n";
                return kExitError;
            }
            const CheckOutcome r = check_circuits(a, b, opts);
            if (opts.json) print_json(out, opts, r);
            else if (!opts.quiet) print_text(out, opts, r);
            return r.equiv.equivalent ? kExitEquivalent : kExitNotEquivalent;
        } catch (const ParseError& e) {
            err << e.what() << "\n";  // already "file:line:column: message"
        } catch (const UsageError& e) {
            err << "qcanopy check: " << e.what() << "\n";
        } catch (const std::exception& e) {
            err << "qcanopy check: " << e.what() << "\n";
        }
        return kExitError;
    }

    int check_main(int argc, const char* const* argv) {
        CheckOptions opts;
        try {
            opts = parse_check_args(argc, argv);
        } catch (const UsageError& e) {
            std::cerr << "qcanopy check: " << e.what() << "\nTry 'qcanopy check --help'.\n";
            return kExitError;
        }
        if (opts.help) { std::cout << check_help; return kExitEquivalent; }
        return run_check(opts, std::cout, std::cerr);
    }

} // namespace qcanopy::cli
