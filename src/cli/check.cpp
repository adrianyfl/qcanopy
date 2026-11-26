#include "check.hpp"

#include "commands.hpp"
#include "qcanopy/parser.hpp"
#include "qcanopy/stimulus.hpp"

#include <cmath>
#include <complex>
#include <exception>
#include <iomanip>
#include <iostream>
#include <set>

namespace qcanopy::cli {

    const char* const check_help =
        "qcanopy check A B [options]\n"
        "Check whether two circuits are equal up to a global phase. A and B are OpenQASM 3.0\n"
        "files, or JSON IR files if their names end in .json. 'qcanopy A B' is short for this.\n"
        "\n"
        "  --stimulus LIST   Comma-separated input states to simulate (default haar,product):\n"
        "                      haar       Haar-random states (--count of them)\n"
        "                      product    random product states (--count of them)\n"
        "                      basis      random computational basis states (--count of them)\n"
        "                      all-basis  every basis state; passing them is a proof (at most 12 qubits)\n"
        "  --count N         States per random kind (default 8).\n"
        "  --seed UINT64     Seed for the random kinds (default 759).\n"
        "  --tol X           Tolerance on the residual norm and phase (default 1e-9).\n"
        "  --no-rewrite      Skip the rewriter and simulate the full circuits.\n"
        "  --overlap         Let the rewriter try overlapping rule substitutions (slower).\n"
        "  --json            Print the result as one JSON object.\n"
        "  --quiet           Print nothing; report the verdict by exit status only.\n"
        "  --help            Show this help.\n"
        "\n"
        "The rewriter (rewriting/run.py pair) first removes common boundaries, inverse pairs and\n"
        "library identities; only what it leaves is simulated. If it cannot run, the full circuits\n"
        "are simulated instead. Environment: QCANOPY_PYTHON (interpreter), QCANOPY_REWRITER\n"
        "(directory holding run.py).\n"
        "\n"
        "Exit status: 0 equivalent, 1 not equivalent, 2 usage, parse or I/O error.\n"
        "A non-equivalent verdict is definitive. An equivalent verdict from random stimuli means\n"
        "no stimulus told the circuits apart; from all-basis or the rewriter it is a proof.\n";

    namespace {
        const std::set<std::string> kStimulusKinds{"haar", "product", "basis", "all-basis"};

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

        bool has_all_basis(const CheckOptions& opts) {
            for (const auto& s : opts.stimuli) if (s == "all-basis") return true;
            return false;
        }

        const char* rewrite_name(RewriteStatus s) {
            switch (s) {
                case RewriteStatus::Proved: return "proved";
                case RewriteStatus::Inconclusive: return "inconclusive";
                case RewriteStatus::Failed: return "failed";
                case RewriteStatus::Skipped: break;
            }
            return "skipped";
        }

        void print_text(std::ostream& out, const CheckOptions& opts, const CheckOutcome& r) {
            const EquivResult& e = r.equiv;
            const bool on_residual = r.rewrite == RewriteStatus::Inconclusive;
            out << std::setprecision(6);
            if (e.equivalent) {
                double phase = std::arg(e.global_phase);
                if (std::abs(phase) < 1e-12) phase = 0.0;  // hide rounding noise like -2.9e-17
                out << "equivalent up to global phase e^{i*" << phase << "}";
                if (r.rewrite == RewriteStatus::Proved) {
                    out << " (proved by the rewriter in " << r.proof_steps << " steps)\n";
                } else {
                    out << " (" << e.stimuli_checked << " stimuli" << (on_residual ? " on the rewriter's residual" : "")
                        << ", max deviation " << std::setprecision(2) << e.deviation << "; "
                        << (has_all_basis(opts) ? "all basis states: proof" : "not a proof") << ")\n";
                }
            } else if (r.rewriter_contradicted) {
                out << "NOT equivalent: simulating the full circuits disagrees with the rewriter "
                    << (r.rewrite == RewriteStatus::Proved ? "proof" : "residual check")
                    << " (deviation " << std::setprecision(2) << e.deviation << "); this is a rewriter bug\n";
            } else {
                out << "NOT equivalent: stimulus #" << e.failing_index << " (" << kind_name(e.failing_kind) << ")"
                    << (on_residual ? " on the rewriter's residual" : "") << ": "
                    << (e.mismatch == Mismatch::State ? "output is not a phase multiple of the input"
                                                      : "output phase differs from stimulus #0's")
                    << " (deviation " << std::setprecision(2) << e.deviation << ")\n";
            }
            if (r.rewrite == RewriteStatus::Proved || on_residual)
                out << "  gates: " << r.gates_a << " vs " << r.gates_b << ", left after rewriting " << r.residual_a
                    << " vs " << r.residual_b << "\n";
        }

        void print_json(std::ostream& out, const CheckOptions& opts, const CheckOutcome& r) {
            const EquivResult& e = r.equiv;
            out << std::setprecision(17);
            out << "{\"a\":" << json_string(opts.file_a) << ",\"b\":" << json_string(opts.file_b)
                << ",\"qubits\":" << r.qubits << ",\"equivalent\":" << (e.equivalent ? "true" : "false")
                << ",\"rewriter\":\"" << rewrite_name(r.rewrite) << "\"";
            if (r.rewrite == RewriteStatus::Failed) out << ",\"rewriter_error\":" << json_string(r.rewrite_error);
            out << ",\"proof_steps\":" << r.proof_steps
                << ",\"rewriter_contradicted\":" << (r.rewriter_contradicted ? "true" : "false")
                << ",\"phase_measured\":" << (r.phase_measured ? "true" : "false")
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
        const Args args(argc, argv, {"--no-rewrite", "--overlap", "--json", "--quiet"},
                        {"--stimulus", "--count", "--seed", "--tol"});
        CheckOptions opts;
        if (args.help()) { opts.help = true; return opts; }
        args.expect_positional(2, 2, "two circuit files");
        opts.file_a = args.positional()[0];
        opts.file_b = args.positional()[1];
        opts.rewrite = !args.has("--no-rewrite");
        opts.rewriter.overlap = args.has("--overlap");
        opts.json = args.has("--json");
        opts.quiet = args.has("--quiet");
        if (auto v = args.get("--stimulus")) opts.stimuli = parse_stimuli(*v);
        if (auto v = args.get("--count")) opts.count = static_cast<std::size_t>(parse_uint("--count", *v));
        if (auto v = args.get("--seed")) opts.seed = parse_uint("--seed", *v);
        if (auto v = args.get("--tol")) opts.tol = parse_real("--tol", *v);

        if (opts.json && opts.quiet) throw UsageError("--json and --quiet cannot be combined");
        if (opts.rewriter.overlap && !opts.rewrite) throw UsageError("--overlap needs the rewriter; drop --no-rewrite");
        if (opts.count == 0) throw UsageError("--count must be at least 1");
        if (!(opts.tol > 0.0)) throw UsageError("--tol must be positive");
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

        const EquivalenceChecker checker(opts.tol);

        // Step 2: let the prover shrink the pair to residuals with A ~ B iff res_A ~ res_B.
        Circuit a_mid = a, b_mid = b;
        if (opts.rewrite) {
            try {
                PairProof proof = prove_pair(a, b, opts.rewriter);
                r.rewrite = proof.proved ? RewriteStatus::Proved : RewriteStatus::Inconclusive;
                r.proof_steps = proof.steps;
                a_mid = std::move(proof.left);
                b_mid = std::move(proof.right);
            } catch (const RewriterError& e) {
                r.rewrite = RewriteStatus::Failed;
                r.rewrite_error = e.what();
            }
        }
        r.residual_a = a_mid.size();
        r.residual_b = b_mid.size();

        // Step 3: simulate whatever is left.
        if (r.rewrite == RewriteStatus::Proved) {
            r.equiv.equivalent = true;
        } else {
            const Stimulus stim = build_stimulus(a.n, opts);
            r.stimuli = stim.size();
            r.equiv = checker.check(a_mid, b_mid, stim);
        }

        // Step 4: residual phases are meaningless, so measure the A-B phase on one state of
        // the full pair. This also catches a wrong proof (a non-equivalent pair fails a
        // Haar-random state with probability 1).
        const bool rewritten = r.rewrite == RewriteStatus::Proved || r.rewrite == RewriteStatus::Inconclusive;
        if (rewritten && r.equiv.equivalent) {
            const EquivResult full = checker.check(a, b, haar_random_states(a.n, 1, opts.seed));
            r.phase_measured = true;
            if (full.equivalent) {
                r.equiv.global_phase = full.global_phase;
            } else {
                r.rewriter_contradicted = true;
                r.equiv = full;
            }
        }
        return r;
    }

    int run_check(const CheckOptions& opts, std::ostream& out, std::ostream& err) {
        try {
            const Circuit a = load_circuit(opts.file_a);
            const Circuit b = load_circuit(opts.file_b);
            if (a.n != b.n) {
                err << "qcanopy check: " << opts.file_a << " has " << a.n << " qubits but " << opts.file_b
                    << " has " << b.n << "\n";
                return kExitError;
            }
            const CheckOutcome r = check_circuits(a, b, opts);
            if (r.rewrite == RewriteStatus::Failed && !opts.quiet)
                err << "qcanopy check: warning: rewriter not used (" << r.rewrite_error
                    << "); simulating the full circuits\n";
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
        if (opts.help) { std::cout << check_help; return kExitOk; }
        return run_check(opts, std::cout, std::cerr);
    }

} // namespace qcanopy::cli
