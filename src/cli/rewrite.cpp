// `qcanopy rewrite`: run one rewriting stage on its own.
//
//   normalize FILE    C++ normalize() (normalize.hpp), prints the normalized circuit
//   peel A B          C++ normalize() + peel() (peel.hpp), prints both residuals
//   prove A B         Python pair prover (rewriter.hpp), prints its verdict
//   classify FILE     Python boundary classifier (rewriting/run.py classify)

#include "commands.hpp"
#include "common.hpp"
#include "rewriter.hpp"

#include "json.hpp"  // src/generator: reading classify output
#include "qcanopy/dag.hpp"
#include "qcanopy/normalize.hpp"
#include "qcanopy/peel.hpp"

#include <iostream>
#include <sstream>
#include <string>

namespace qcanopy::cli {

    namespace {
        const char* const help =
            "qcanopy rewrite <subcommand> ...\n"
            "Run one rewriting stage on its own. Files are OpenQASM 3.0, or JSON IR v1 if *.json.\n"
            "\n"
            "  normalize FILE [--emit qasm|ir]\n"
            "      Canonicalize gate spellings, fuse rotations and cancel neighbouring inverses\n"
            "      (exact, including global phase). Prints the result (default qasm).\n"
            "  peel A B [--emit qasm|ir] [--no-normalize]\n"
            "      Normalize both, then strip their common prefix and suffix. Prints what is left\n"
            "      of each; peeled phases are added to A's global phase.\n"
            "  prove A B [--overlap] [--json]\n"
            "      Run the Python pair prover (rewriting/run.py pair). Exit 0 proved, 1 inconclusive.\n"
            "      Inconclusive does not mean the circuits differ; use 'qcanopy check' for a verdict.\n"
            "  classify FILE [--overlap] [--json]\n"
            "      Run the Python boundary classifier (rewriting/run.py classify). Exit 0 if the\n"
            "      search completed within its bounds, 1 if it stopped early.\n"
            "\n"
            "--json prints the Python tool's full JSON output. The Python stages use\n"
            "QCANOPY_PYTHON (interpreter) and QCANOPY_REWRITER (directory holding run.py).\n"
            "Exit status 2 means a usage, parse, I/O or rewriter error.\n";

        Format emit_format(const Args& args) {
            const auto v = args.get("--emit");
            return v ? parse_format("--emit", *v) : Format::Qasm;
        }

        void same_width(const Circuit& a, const Circuit& b) {
            if (a.n != b.n)
                throw std::invalid_argument("circuits have different qubit counts (" + std::to_string(a.n) +
                                            " and " + std::to_string(b.n) + ")");
        }

        int normalize_cmd(int argc, const char* const* argv) {
            const Args args(argc, argv, {}, {"--emit"});
            if (args.help()) { std::cout << help; return kExitOk; }
            args.expect_positional(1, 1, "one circuit file");
            const Format format = emit_format(args);
            CircuitDAG dag(load_circuit(args.positional()[0]));
            normalize(dag);
            write_circuit(std::cout, dag.to_circuit(), format);
            return kExitOk;
        }

        int peel_cmd(int argc, const char* const* argv) {
            const Args args(argc, argv, {"--no-normalize"}, {"--emit"});
            if (args.help()) { std::cout << help; return kExitOk; }
            args.expect_positional(2, 2, "two circuit files");
            const Format format = emit_format(args);
            const std::string& fa = args.positional()[0];
            const std::string& fb = args.positional()[1];
            const Circuit ca = load_circuit(fa), cb = load_circuit(fb);
            same_width(ca, cb);

            CircuitDAG a(ca), b(cb);
            if (!args.has("--no-normalize")) { normalize(a); normalize(b); }
            std::size_t pairs = 0;
            for (std::size_t k; (k = peel(/*R=*/b, /*C=*/a)) > 0;) pairs += k;  // phases go to A
            const Circuit ra = a.to_circuit(), rb = b.to_circuit();

            if (format == Format::Ir) {
                auto ir = [](const Circuit& c) {
                    std::ostringstream s;
                    write_circuit(s, c, Format::Ir);
                    std::string text = s.str();
                    while (!text.empty() && text.back() == '\n') text.pop_back();
                    return text;
                };
                std::cout << "{\"peeled_pairs\":" << pairs << ",\"a\":" << ir(ra) << ",\"b\":" << ir(rb) << "}\n";
            } else {
                std::cout << "// qcanopy rewrite peel: removed " << pairs << " gate pairs; " << fa << " has "
                          << ra.size() << " gates left, " << fb << " has " << rb.size() << "\n"
                          << "// ---- residual of " << fa << " ----\n";
                write_circuit(std::cout, ra, Format::Qasm);
                std::cout << "// ---- residual of " << fb << " ----\n";
                write_circuit(std::cout, rb, Format::Qasm);
            }
            return kExitOk;
        }

        int prove_cmd(int argc, const char* const* argv) {
            const Args args(argc, argv, {"--overlap", "--json"}, {});
            if (args.help()) { std::cout << help; return kExitOk; }
            args.expect_positional(2, 2, "two circuit files");
            const Circuit a = load_circuit(args.positional()[0]), b = load_circuit(args.positional()[1]);
            same_width(a, b);
            RewriterConfig cfg = default_rewriter_config();
            cfg.overlap = args.has("--overlap");

            const RewriterRun run = run_rewriter(cfg, "pair", {{"--left", &a}, {"--right", &b}});
            std::istringstream in(run.output);
            const PairProof p = parse_pair_proof(in, a.n, /*convert_residuals=*/false);
            if (args.has("--json")) {
                std::cout << run.output;
            } else {
                std::cout << (p.proved ? "proved: equivalent up to global phase" : "inconclusive") << " ("
                          << p.steps << " proof steps; " << a.size() << " vs " << b.size() << " gates, "
                          << p.residual_gates_left << " vs " << p.residual_gates_right << " left)\n";
                if (!p.step_kinds.empty()) {
                    std::cout << "  steps:";
                    for (const auto& k : p.step_kinds) std::cout << ' ' << k;
                    std::cout << "\n";
                }
            }
            return p.proved ? kExitOk : 1;
        }

        int classify_cmd(int argc, const char* const* argv) {
            const Args args(argc, argv, {"--overlap", "--json"}, {});
            if (args.help()) { std::cout << help; return kExitOk; }
            args.expect_positional(1, 1, "one circuit file");
            const Circuit c = load_circuit(args.positional()[0]);
            RewriterConfig cfg = default_rewriter_config();
            cfg.overlap = args.has("--overlap");

            const RewriterRun run = run_rewriter(cfg, "classify", {{"--dag", &c}});
            if (args.has("--json")) {
                std::cout << run.output;
            } else {
                std::istringstream in(run.output);
                const auto root = testgen::json::Parser(in).parse();
                const auto& cls = root.field("classification");
                const bool complete = cls.field("complete_for_bounds").number != 0;
                std::cout << c.size() << " gates: " << cls.field("matches").items().size() << " matches, "
                          << cls.field("regions").items().size() << " regions, "
                          << cls.field("composition_groups").items().size() << " composition groups, "
                          << cls.field("unmatched_gate_ids").items().size() << " gates unmatched\n"
                          << "  search " << (complete ? "complete" : "stopped early") << " ("
                          << cls.field("stop_reason").str() << ")\n";
                if (cfg.overlap && cls.object.count("rewrite_candidates"))
                    std::cout << "  overlap: " << cls.field("rewrite_candidates").items().size()
                              << " rewrite candidates\n";
            }
            return run.status == 0 ? kExitOk : 1;
        }
    } // namespace

    int rewrite_main(int argc, const char* const* argv) {
        const std::string sub = argc >= 2 ? argv[1] : "";
        if (sub == "--help" || sub == "-h") { std::cout << help; return kExitOk; }
        if (sub.empty()) { std::cerr << help; return kExitError; }

        int (*cmd)(int, const char* const*) = nullptr;
        if (sub == "normalize") cmd = normalize_cmd;
        else if (sub == "peel") cmd = peel_cmd;
        else if (sub == "prove") cmd = prove_cmd;
        else if (sub == "classify") cmd = classify_cmd;
        else {
            std::cerr << "qcanopy rewrite: unknown subcommand '" << sub
                      << "' (expected normalize, peel, prove or classify)\nTry 'qcanopy rewrite --help'.\n";
            return kExitError;
        }
        return guarded("rewrite " + sub, [&] { return cmd(argc - 1, argv + 1); });
    }

} // namespace qcanopy::cli
