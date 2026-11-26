// qcanopy — command-line front end. Each subcommand lives in its own file under src/cli/;
// this file only dispatches on the first argument.

#include "check.hpp"

#include <cstdio>
#include <string>

namespace {
    const char* const usage =
        "usage: qcanopy <command> [args]\n"
        "\n"
        "Commands:\n"
        "  check A.qasm B.qasm [options]   Check two circuits for equivalence up to global phase\n"
        "\n"
        "Run 'qcanopy <command> --help' for a command's options.\n";
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fputs(usage, stderr); return qcanopy::cli::kExitError; }

    const std::string cmd = argv[1];
    if (cmd == "check") return qcanopy::cli::check_main(argc - 1, argv + 1);
    if (cmd == "--help" || cmd == "-h" || cmd == "help") { std::fputs(usage, stdout); return 0; }

    std::fprintf(stderr, "qcanopy: unknown command '%s'\n\n%s", cmd.c_str(), usage);
    return qcanopy::cli::kExitError;
}
