// src/cli/commands.hpp — entry points of the `qcanopy` subcommands.
//
// Each takes the arguments starting at the command name (argv[0] is skipped) and returns
// an exit status. `check` is the main command: parse, rewrite and simulate in one step. The
// others expose one stage each:
//
//   check A B            check.cpp     equivalence check (the whole pipeline)
//   parse FILE...        parse.cpp     summarize, or convert between OpenQASM and JSON IR
//   rewrite SUB ...      rewrite.cpp   normalize | peel (C++), prove | classify (Python)
//   sim FILE             sim.cpp       run a circuit and print the output state
//   generate ...         generate.cpp  labeled test-case pairs (dataset_io.hpp)

#pragma once

namespace qcanopy::cli {

    int check_main(int argc, const char* const* argv);
    int parse_main(int argc, const char* const* argv);
    int rewrite_main(int argc, const char* const* argv);
    int sim_main(int argc, const char* const* argv);
    int generate_main(int argc, const char* const* argv);

} // namespace qcanopy::cli
