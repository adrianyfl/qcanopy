// src/cli/common.hpp — pieces shared by every `qcanopy` subcommand.
//
// Exit status, unless a command's help says otherwise: 0 success, 2 usage, parse or I/O error.
// Circuit files are OpenQASM 3.0 (parser.hpp), or JSON IR v1 (dataset_io.hpp) when the name
// ends in ".json".

#pragma once

#include "qcanopy/ir.hpp"

#include <cstdint>
#include <functional>
#include <initializer_list>
#include <iosfwd>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace qcanopy::cli {

    inline constexpr int kExitOk = 0;
    inline constexpr int kExitEquivalent = 0;
    inline constexpr int kExitNotEquivalent = 1;
    inline constexpr int kExitError = 2;

    // Bad command line. what() is the message without the "qcanopy <command>:" prefix.
    class UsageError : public std::invalid_argument {
    public:
        using std::invalid_argument::invalid_argument;
    };

    // The arguments of one subcommand, split into flags, "--key value" options and positionals.
    // "--help" is always accepted and stops parsing (so it wins over later bad arguments).
    // Unknown, duplicate and value-less options throw UsageError.
    class Args {
    public:
        // argv[0] is the command name and is skipped.
        Args(int argc, const char* const* argv, std::initializer_list<const char*> flags,
             std::initializer_list<const char*> valued);

        bool help() const { return help_; }
        const std::vector<std::string>& positional() const { return positional_; }
        bool has(const std::string& key) const { return flags_.count(key) || values_.count(key); }
        std::optional<std::string> get(const std::string& key) const;

        // Throws UsageError unless there are between min and max positionals.
        void expect_positional(std::size_t min, std::size_t max, const char* what) const;

    private:
        bool help_ = false;
        std::vector<std::string> positional_;
        std::set<std::string> flags_;
        std::map<std::string, std::string> values_;
    };

    // Strict number parsing for option values; `key` names the option in the error.
    uint64_t parse_uint(const std::string& key, const std::string& s);
    double parse_real(const std::string& key, const std::string& s);

    enum class Format : uint8_t { Qasm, Ir };
    Format parse_format(const std::string& key, const std::string& s);   // "qasm" | "ir"

    // Reads a circuit (see the header comment). Throws ParseError ("file:line:column: ...")
    // for bad QASM and std::runtime_error ("file: ...") otherwise.
    Circuit load_circuit(const std::string& path);

    // Writes c as OpenQASM 3.0 or JSON IR v1 (one object and a newline).
    void write_circuit(std::ostream& out, const Circuit& c, Format format);

    // s as a JSON string literal, quotes included.
    std::string json_string(const std::string& s);

    // Runs body, turning exceptions into a message on stderr and kExitError. `command` is the
    // name used in messages, e.g. "parse" or "rewrite peel".
    int guarded(const std::string& command, const std::function<int()>& body);

} // namespace qcanopy::cli
