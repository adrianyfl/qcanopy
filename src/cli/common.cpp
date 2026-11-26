#include "common.hpp"

#include "json.hpp"                // src/generator: json::quote
#include "qcanopy/dataset_io.hpp"  // read_ir, write_ir, write_qasm
#include "qcanopy/parser.hpp"

#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <locale>
#include <new>
#include <sstream>

namespace qcanopy::cli {

    Args::Args(int argc, const char* const* argv, std::initializer_list<const char*> flags,
               std::initializer_list<const char*> valued) {
        const std::set<std::string> flag_set(flags.begin(), flags.end());
        const std::set<std::string> valued_set(valued.begin(), valued.end());
        for (int i = 1; i < argc; ++i) {
            const std::string key = argv[i];
            if (key == "--help" || key == "-h") { help_ = true; return; }
            if (key.size() < 2 || key.rfind("--", 0) != 0) { positional_.push_back(key); continue; }
            if (has(key)) throw UsageError("duplicate option " + key);
            if (flag_set.count(key)) { flags_.insert(key); continue; }
            if (!valued_set.count(key)) throw UsageError("unknown option " + key);
            if (i + 1 == argc) throw UsageError("missing value for " + key);
            values_[key] = argv[++i];
        }
    }

    std::optional<std::string> Args::get(const std::string& key) const {
        const auto it = values_.find(key);
        if (it == values_.end()) return std::nullopt;
        return it->second;
    }

    void Args::expect_positional(std::size_t min, std::size_t max, const char* what) const {
        const std::size_t n = positional_.size();
        if (n >= min && n <= max) return;
        throw UsageError(std::string("expected ") + what + ", got " + std::to_string(n) + " argument" +
                         (n == 1 ? "" : "s"));
    }

    uint64_t parse_uint(const std::string& key, const std::string& s) {
        uint64_t out = 0;
        const auto r = std::from_chars(s.data(), s.data() + s.size(), out);
        if (s.empty() || r.ec != std::errc{} || r.ptr != s.data() + s.size())
            throw UsageError(key + " expects a non-negative integer, got '" + s + "'");
        return out;
    }

    double parse_real(const std::string& key, const std::string& s) {
        std::istringstream in(s);
        in.imbue(std::locale::classic());
        double out = 0.0;
        in >> std::noskipws >> out;
        if (!in || !in.eof() || !std::isfinite(out))
            throw UsageError(key + " expects a finite number, got '" + s + "'");
        return out;
    }

    Format parse_format(const std::string& key, const std::string& s) {
        if (s == "qasm") return Format::Qasm;
        if (s == "ir") return Format::Ir;
        throw UsageError(key + " expects qasm or ir, got '" + s + "'");
    }

    Circuit load_circuit(const std::string& path) {
        if (std::filesystem::path(path).extension() != ".json") return parse_qasm_file(path);
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error(path + ": cannot read file");
        try {
            return testgen::read_ir(in);
        } catch (const std::exception& e) {
            throw std::runtime_error(path + ": " + e.what());
        }
    }

    void write_circuit(std::ostream& out, const Circuit& c, Format format) {
        if (format == Format::Ir) testgen::write_ir(out, c);
        else testgen::write_qasm(out, c);
    }

    std::string json_string(const std::string& s) { return testgen::json::quote(s); }

    int guarded(const std::string& command, const std::function<int()>& body) {
        try {
            return body();
        } catch (const UsageError& e) {
            std::cerr << "qcanopy " << command << ": " << e.what() << "\nTry 'qcanopy " << command << " --help'.\n";
        } catch (const ParseError& e) {
            std::cerr << e.what() << "\n";  // already "file:line:column: message"
        } catch (const std::bad_alloc&) {
            std::cerr << "qcanopy " << command << ": out of memory\n";
        } catch (const std::exception& e) {
            std::cerr << "qcanopy " << command << ": " << e.what() << "\n";
        }
        return kExitError;
    }

} // namespace qcanopy::cli
