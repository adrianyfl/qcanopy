// src/cli/check.hpp — the `qcanopy check` subcommand.
//
//   qcanopy check A.qasm B.qasm [options]
//
// 1. Both files are loaded (OpenQASM, or JSON IR if *.json; common.hpp) and must have the
//    same qubit count.
// 2. Unless --no-rewrite is given, the pair goes to the Python pair prover (rewriter.hpp),
//    which removes common boundaries, gate/inverse pairs and library identities. If it proves
//    the pair, no stimuli are simulated. If it is inconclusive, only the residual it leaves
//    is simulated: A ~ B iff residual_A ~ residual_B. If it fails (no Python, prover error,
//    a residual gate the IR cannot hold), a warning is printed and step 3 checks A and B.
// 3. The (residual) circuits are compared by EquivalenceChecker (equiv.hpp) on the stimulus
//    batch selected by --stimulus / --count / --seed.
// 4. The prover does not track global phase. So whenever it was used and the verdict is
//    equivalent, A and B themselves are simulated on one Haar-random state to measure the
//    phase. That run also guards against a wrong proof: if it disagrees, the verdict is
//    NOT equivalent and rewriter_contradicted is set.
//
// Exit status follows cmp(1): 0 equivalent, 1 not equivalent, 2 usage, parse or I/O error.
//
// The logic lives here rather than in main.cpp so tests can call parse_check_args() and
// check_circuits() directly.

#pragma once

#include "common.hpp"
#include "qcanopy/equiv.hpp"
#include "qcanopy/ir.hpp"
#include "rewriter.hpp"

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <vector>

namespace qcanopy::cli {

    // all-basis needs 2^n states of 2^n amplitudes (n = 12 -> 256 MiB), so it is refused above this.
    inline constexpr int kMaxAllBasisQubits = 12;

    struct CheckOptions {
        std::string file_a, file_b;
        std::vector<std::string> stimuli{"haar", "product"};  // haar | product | basis | all-basis
        std::size_t count = 8;                                 // states per random kind
        uint64_t seed = 759;
        double tol = 1e-9;
        bool rewrite = true;                                   // run the pair prover first
        RewriterConfig rewriter = default_rewriter_config();   // interpreter, run.py dir, --overlap
        bool json = false;
        bool quiet = false;
        bool help = false;                                     // --help was given; nothing else is valid
    };

    enum class RewriteStatus : uint8_t {
        Skipped,        // --no-rewrite
        Proved,         // the prover removed everything
        Inconclusive,   // a residual was left and simulated
        Failed,         // the prover could not be used; A and B were simulated instead
    };

    struct CheckOutcome {
        EquivResult equiv;                 // verdict (step 3, with the phase from step 4)
        RewriteStatus rewrite = RewriteStatus::Skipped;
        std::string rewrite_error;         // why, if rewrite == Failed
        std::size_t proof_steps = 0;       // steps the prover recorded
        bool phase_measured = false;       // step 4 ran
        bool rewriter_contradicted = false;  // step 4 found the pair NOT equivalent after all
        int qubits = 0;
        std::size_t gates_a = 0, gates_b = 0;          // as parsed
        std::size_t residual_a = 0, residual_b = 0;    // simulated in step 3 (= gates_* unless rewritten)
        std::size_t stimuli = 0;                       // batch size handed to the checker in step 3
    };

    extern const char* const check_help;

    // Parses the arguments after "check" (argv[0] is "check"). Throws UsageError.
    CheckOptions parse_check_args(int argc, const char* const* argv);

    // Builds the stimulus batch for n qubits. Throws UsageError for an unknown kind or for
    // all-basis above kMaxAllBasisQubits.
    Stimulus build_stimulus(int n, const CheckOptions& opts);

    // Runs steps 2 to 4 above. Throws std::invalid_argument if a.n != b.n.
    CheckOutcome check_circuits(const Circuit& a, const Circuit& b, const CheckOptions& opts);

    // Loads opts.file_a / opts.file_b, checks them, and writes the report to `out`
    // (errors to `err`). Returns an exit status.
    int run_check(const CheckOptions& opts, std::ostream& out, std::ostream& err);

    // check_main() is declared in commands.hpp.

} // namespace qcanopy::cli
