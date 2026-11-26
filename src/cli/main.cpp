// qcanopy — command-line front end. `check` is the main command (the whole pipeline); the
// others expose one stage each (commands.hpp). This file only dispatches on the first argument.

#include "commands.hpp"
#include "common.hpp"

#include <cstdio>
#include <filesystem>
#include <string>

#ifndef QCANOPY_VERSION
#define QCANOPY_VERSION "unknown"
#endif

namespace {
    const char* const usage =
        "usage: qcanopy <command> [args]\n"
        "       qcanopy A B [options]          short for 'qcanopy check A B'\n"
        "\n"
        "Main command:\n"
        "  check A B       Check two circuits for equivalence up to global phase:\n"
        "                  parse, rewrite, then simulate whatever the rewriter leaves\n"
        "\n"
        "Single stages:\n"
        "  parse FILE...   Summarize circuits, or convert between OpenQASM and JSON IR\n"
        "  rewrite SUB     normalize | peel (C++ passes), prove | classify (Python rewriter)\n"
        "  sim FILE        Simulate a circuit and print its output state\n"
        "  generate ...    Generate labeled equivalent / non-equivalent test-case pairs\n"
        "\n"
        "Circuit files are OpenQASM 3.0, or JSON IR v1 if their names end in .json.\n"
        "Run 'qcanopy <command> --help' for a command's options; 'qcanopy --version' for the version.\n";

    using Command = int (*)(int, const char* const*);

    Command find(const std::string& name) {
        using namespace qcanopy::cli;
        if (name == "check") return check_main;
        if (name == "parse") return parse_main;
        if (name == "rewrite") return rewrite_main;
        if (name == "sim") return sim_main;
        if (name == "generate") return generate_main;
        return nullptr;
    }
}

int main(int argc, char** argv) {
    using qcanopy::cli::kExitError;
    if (argc < 2) { std::fputs(usage, stderr); return kExitError; }

    const std::string first = argv[1];
    if (first == "--help" || first == "-h" || first == "help") { std::fputs(usage, stdout); return 0; }
    if (first == "--version") { std::puts("qcanopy " QCANOPY_VERSION); return 0; }
    if (const Command cmd = find(first)) return cmd(argc - 1, argv + 1);

    // `qcanopy A B ...`: a file where the command should be means check. argv[0] stands in
    // for the command name, which check skips.
    std::error_code ec;
    if (std::filesystem::is_regular_file(first, ec)) return qcanopy::cli::check_main(argc, argv);

    std::fprintf(stderr, "qcanopy: unknown command '%s'\n\n%s", first.c_str(), usage);
    return kExitError;
}
