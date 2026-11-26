// src/cli/rewriter.hpp — run the Python pair prover (rewriting/, see rewriting/README.md).
//
// The rewriting module has no in-process C++ bindings, so this runs
//
//     <python> <dir>/run.py pair --left A.json --right B.json --input-format qcanopy-ir -o out.json
//
// as a child process (no shell), exchanging circuits as JSON IR v1 (dataset_io.hpp, write_ir).
// run_rewriter() runs any run.py command the same way (`qcanopy rewrite classify` uses it).
//
// 1. What comes back. The prover removes common boundary regions, gate/inverse pairs and
//    library identities from both circuits. Its verdict is "equivalent_up_to_global_phase"
//    (nothing left) or "inconclusive", in which case the residual DAGs are returned.
//    Every step replaces a region by one equal to it UP TO PHASE, so A ~ B iff residual_A ~
//    residual_B, but the phase between A and B is NOT tracked: residual global phases are
//    meaningless and are set to 0 here.
//
// 2. Coordinates. If the prover compared SWAP-free cores (layout.cores_used), both residuals
//    are on core wires with the same final permutation, so the equivalence above still holds.
//
// 3. Failure. Any problem (no interpreter, prover error, unexpected output, or a residual
//    gate the IR cannot represent, e.g. a library gate introduced with --overlap) throws
//    RewriterError. Callers should then fall back to checking the original circuits.
//
// Only implemented on POSIX systems; elsewhere prove_pair() always throws RewriterError.

#pragma once

#include "qcanopy/ir.hpp"

#include <cstddef>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qcanopy::cli {

    class RewriterError : public std::runtime_error {
    public:
        using std::runtime_error::runtime_error;
    };

    struct RewriterConfig {
        std::string python;       // interpreter (looked up on PATH if it has no '/')
        std::string dir;          // directory holding run.py
        bool overlap = false;     // pass --enable-overlap
    };

    // $QCANOPY_PYTHON and $QCANOPY_REWRITER if set, else the paths found at configure time.
    RewriterConfig default_rewriter_config();

    struct PairProof {
        bool proved = false;          // verdict was equivalent_up_to_global_phase
        bool cores_used = false;      // residuals are on SWAP-free core wires (point 2)
        std::size_t steps = 0;        // proof steps the prover recorded
        std::vector<std::string> step_kinds;   // "structural", "inverse_identity", ...
        std::size_t residual_gates_left = 0, residual_gates_right = 0;
        Circuit left, right;          // residuals of A and B, global_phase 0 (point 1);
                                      // empty circuits unless convert_residuals was set
    };

    // Runs the prover on A and B (same qubit count). Throws RewriterError (point 3). With
    // convert_residuals false, residual gates are only counted, so any gate is accepted.
    PairProof prove_pair(const Circuit& a, const Circuit& b, const RewriterConfig& cfg,
                         bool convert_residuals = true);

    // Reads the prover's qrewrite-pair/1 output for n-qubit circuits. Throws RewriterError.
    // Exposed for tests; prove_pair() calls it on the child's output.
    PairProof parse_pair_proof(std::istream& json, int n, bool convert_residuals = true);

    // What run.py printed: its exit status (0 completed/proved, 2 partial/inconclusive) and
    // its JSON output.
    struct RewriterRun {
        int status = 0;
        std::string output;
    };

    // Runs `run.py <command> <option> <file>... --input-format qcanopy-ir -o <out>`, writing
    // each input circuit to a temp file passed after its option (e.g. {"--dag", &c}). Throws
    // RewriterError on any other exit status or if the rewriter cannot be started.
    RewriterRun run_rewriter(const RewriterConfig& cfg, const std::string& command,
                             const std::vector<std::pair<std::string, const Circuit*>>& inputs);

} // namespace qcanopy::cli
