#pragma once
#include "qcanopy/generator.hpp"
#include <filesystem>
#include <iosfwd>

namespace qcanopy::testgen {
namespace fs = std::filesystem;
struct OutputOptions {
    std::string format;
    fs::path destination;
};
void write_ir(std::ostream&, const Circuit&);
Circuit read_ir(std::istream&);
void write_qasm(std::ostream&, const Circuit&);
void validate_output_options(const OutputOptions&);
// Publishes manifest.jsonl last. Refuses nonempty destinations; never overwrites files.
void write_dataset(const OutputOptions&, const GeneratorConfig&, const std::vector<CircuitPair>&);
} // namespace qcanopy::testgen
