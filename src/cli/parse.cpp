// `qcanopy parse`: read circuit files and summarize them, or convert one between formats.
// The summary output and exit statuses match the older qcanopy_parse tool, which now runs this.

#include "commands.hpp"
#include "common.hpp"

#include <cstdio>
#include <iostream>
#include <string>

namespace qcanopy::cli {

    namespace {
        const char* const help =
            "qcanopy parse FILE... [--emit qasm|ir]\n"
            "Parse each circuit file and print its qubit count, gate count and global phase.\n"
            "FILE is OpenQASM 3.0, or JSON IR v1 if its name ends in .json.\n"
            "\n"
            "  --emit qasm|ir   Instead of the summary, print the single FILE as OpenQASM 3.0 or\n"
            "                   JSON IR v1 (converts between the two formats).\n"
            "  --help           Show this help.\n"
            "\n"
            "A file that cannot be read or parsed is reported on stderr; parse errors read\n"
            "FILE:LINE:COLUMN: message.\n"
            "Exit status: 0 if every file parsed, 1 if any failed, 2 on a usage error.\n";
    }

    int parse_main(int argc, const char* const* argv) {
        return guarded("parse", [&] {
            const Args args(argc, argv, {}, {"--emit"});
            if (args.help()) { std::cout << help; return kExitOk; }
            args.expect_positional(1, SIZE_MAX, "at least one circuit file");

            if (const auto emit = args.get("--emit")) {
                const Format format = parse_format("--emit", *emit);
                args.expect_positional(1, 1, "one circuit file with --emit");
                write_circuit(std::cout, load_circuit(args.positional()[0]), format);
                return kExitOk;
            }

            const auto& files = args.positional();
            int failed = 0;
            for (const auto& file : files) {
                try {
                    const Circuit c = load_circuit(file);
                    std::printf("%s: qubits=%d gates=%zu global_phase=%.17g\n", file.c_str(), c.n, c.size(),
                                c.global_phase);
                } catch (const std::exception& e) {
                    std::fflush(stdout);
                    std::fprintf(stderr, "%s\n", e.what());
                    ++failed;
                }
            }
            if (files.size() > 1)
                std::printf("parsed %zu of %zu files\n", files.size() - static_cast<std::size_t>(failed), files.size());
            return failed ? 1 : kExitOk;
        });
    }

} // namespace qcanopy::cli
