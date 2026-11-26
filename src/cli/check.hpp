// src/cli/check.hpp — the `qcanopy check` subcommand.
//
//   qcanopy check A.qasm B.qasm [options]
//
// 1. Both files are parsed (parser.hpp) and must have the same qubit count.
// 2. Unless --no-rewrite is given, both circuits are normalized and their common prefix and
//    suffix are peeled off (normalize.hpp, peel.hpp). If nothing is left, the pair is proved
//    equivalent without simulation. The peeled phases are carried in the residuals'
//    global_phase, so the reported phase is still the one between A and B.
// 3. Otherwise the residuals are compared by EquivalenceChecker (equiv.hpp) on the stimulus
//    batch selected by --stimulus / --count / --seed.
//
// Exit status follows cmp(1): 0 equivalent, 1 not equivalent, 2 usage, parse or I/O error.
//
// The logic lives here rather than in main.cpp so tests can call parse_check_args() and
// check_circuits() directly.

#pragma once

#include "qcanopy/equiv.hpp"
#include "qcanopy/ir.hpp"

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <vector>

namespace qcanopy::cli {

    // Exit statuses shared by every subcommand.
    inline constexpr int kExitEquivalent = 0;
    inline constexpr int kExitNotEquivalent = 1;
    inline constexpr int kExitError = 2;

    // Bad command line. what() is the message without the "qcanopy check:" prefix.
    class UsageError : public std::invalid_argument {
    public:
        using std::invalid_argument::invalid_argument;
    };

    // all-basis needs 2^n states of 2^n amplitudes (n = 12 -> 256 MiB), so it is refused above this.
    inline constexpr int kMaxAllBasisQubits = 12;

    struct CheckOptions {
        std::string file_a, file_b;
        std::vector<std::string> stimuli{"haar", "product"};  // haar | product | basis | all-basis
        std::size_t count = 8;                                 // states per random kind
        uint64_t seed = 759;
        double tol = 1e-9;
        bool rewrite = true;                                   // normalize + peel pre-pass
        bool json = false;
        bool quiet = false;
        bool help = false;                                     // --help was given; nothing else is valid
    };

    struct CheckOutcome {
        EquivResult equiv;                 // verdict; stimuli_checked == 0 if proved by rewriting
        bool proved_by_rewrite = false;    // the pre-pass left nothing to simulate
        int qubits = 0;
        std::size_t gates_a = 0, gates_b = 0;          // as parsed
        std::size_t residual_a = 0, residual_b = 0;    // after the pre-pass (= gates_* without it)
        std::size_t stimuli = 0;                       // batch size handed to the checker
    };

    extern const char* const check_help;

    // Parses the arguments after "check" (argv[0] is "check"). Throws UsageError.
    CheckOptions parse_check_args(int argc, const char* const* argv);

    // Builds the stimulus batch for n qubits. Throws UsageError for an unknown kind or for
    // all-basis above kMaxAllBasisQubits.
    Stimulus build_stimulus(int n, const CheckOptions& opts);

    // Runs steps 2 and 3 above. Throws std::invalid_argument if a.n != b.n.
    CheckOutcome check_circuits(const Circuit& a, const Circuit& b, const CheckOptions& opts);

    // Loads opts.file_a / opts.file_b, checks them, and writes the report to `out`
    // (errors to `err`). Returns an exit status.
    int run_check(const CheckOptions& opts, std::ostream& out, std::ostream& err);

    // Entry point from main.cpp: parse, run, report. Returns an exit status.
    int check_main(int argc, const char* const* argv);

} // namespace qcanopy::cli
