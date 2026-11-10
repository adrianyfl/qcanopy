#pragma once

#include "qcanopy/ir.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace qcanopy::testgen {
inline constexpr double kPrefixTolerance = 1e-12;
inline constexpr double kVerificationTolerance = 1e-8;
inline constexpr int kDenseQubitCap = 5;
inline constexpr const char* kGeneratorVersion = "1.0.0";

struct GeneratorConfig {
    int qubits = 4;
    std::size_t gates = 100, pairs = 32, prefix_group_size = 1;
    double prefix_fraction = 0.5, equivalent_fraction = 0.5;
    uint64_t seed = 759;
    std::string gate_set = "mixed", angle_mode = "mixed", verify = "small";
};
enum class ExpectedRelation { Equivalent, NonEquivalent };
enum class TransformKind { Identical, InsertInverse, XtoHZH, StoTT, SwapToCX, Commute, GlobalPhase };
const char* transform_name(TransformKind);
const char* relation_name(ExpectedRelation);
struct TransformRecord {
    TransformKind kind = TransformKind::Identical;
    std::size_t location = 0;
    double phase_adjustment = 0;
    std::vector<Gate> removed, inserted;
    std::string proof;
};
struct PrefixMetrics {
    std::size_t gates = 0, a_gates = 0, b_gates = 0;
    double fraction = 0;
};
struct SharingMetrics { PrefixMetrics raw, normalized; };
struct CaseMetadata {
    std::string case_id, prefix_identity, verification;
    std::size_t prefix_group_id = 0, planted_prefix_gates = 0;
    uint64_t master_seed = 0, prefix_seed = 0, suffix_seed = 0, transform_seed = 0;
    TransformRecord transform;
    int negative_qubit = -1;
    SharingMetrics sharing;
};
struct CircuitPair {
    Circuit a, b;
    ExpectedRelation expected = ExpectedRelation::Equivalent;
    CaseMetadata metadata;
};

void validate_config(const GeneratorConfig&);
// Stronger checks than upstream validate: finite params, factory consistency, unitarity.
void validate_circuit(const Circuit&);
Circuit generate_circuit(const GeneratorConfig&, uint64_t seed);
CircuitPair generate_pair(const GeneratorConfig&, std::size_t case_index);
std::vector<CircuitPair> generate_batch(const GeneratorConfig&);
Circuit inverse_circuit(const Circuit&);
TransformRecord apply_equivalent_transform(Circuit&, TransformKind, std::size_t protected_prefix, uint64_t seed);
PrefixMetrics common_prefix(const Circuit&, const Circuit&, double tolerance = kPrefixTolerance);
Circuit normalized_copy(const Circuit&);
SharingMetrics measure_sharing(const Circuit&, const Circuit&);
// Full operator, column-major, including global phase. Throws above kDenseQubitCap.
std::vector<cplx> reference_unitary(const Circuit&);
bool equivalent_up_to_phase(const Circuit&, const Circuit&, double tolerance = kVerificationTolerance);
void verify_pair(CircuitPair&, const std::string& mode = "small");
} // namespace qcanopy::testgen
