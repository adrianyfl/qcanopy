// qcanopy_parse — parses OpenQASM 3.0 files and reports what was read from each.

#include "qcanopy/parser.hpp"

#include <cstdio>
#include <exception>
#include <string>

namespace {
    const char* help =
        "qcanopy_parse FILE...\n"
        "Parse each OpenQASM 3.0 file and print its qubit count, gate count and global phase.\n"
        "A file that cannot be read or parsed is reported on stderr; parse errors read\n"
        "FILE:LINE:COLUMN: message.\n"
        "Exit status: 0 if every file parsed, 1 if any failed, 2 if no file was given.\n";
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fputs(help, stderr); return 2; }
    if (std::string(argv[1]) == "--help") { std::fputs(help, stdout); return 0; }

    const int files = argc - 1;
    int failed = 0;
    for (int i = 1; i < argc; ++i) {
        try {
            const qcanopy::Circuit c = qcanopy::parse_qasm_file(argv[i]);
            std::printf("%s: qubits=%d gates=%zu global_phase=%.17g\n", argv[i], c.n, c.size(), c.global_phase);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "%s\n", e.what());
            ++failed;
        }
    }
    if (files > 1) std::printf("parsed %d of %d files\n", files - failed, files);
    return failed ? 1 : 0;
}
