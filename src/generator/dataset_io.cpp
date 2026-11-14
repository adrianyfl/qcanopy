#include "qcanopy/dataset_io.hpp"
#include "json.hpp"
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

#ifndef QCANOPY_SOURCE_COMMIT
#define QCANOPY_SOURCE_COMMIT "unknown"
#endif

namespace qcanopy::testgen {
namespace {
std::string gate_name(GateKind k) {
    if (k == GateKind::U1Q) return "unitary_1q";
    if (k == GateKind::U2Q) return "unitary_2q";
    return std::string(info(k).name);
}
GateKind gate_kind(const std::string& name) {
    for (int i = 0; i < static_cast<int>(GateKind::Count_); ++i)
        if (gate_name(static_cast<GateKind>(i)) == name) return static_cast<GateKind>(i);
    throw std::invalid_argument("unknown IR gate name: " + name);
}
void gate_json(std::ostream& out, const Gate& g) {
    out << "{\"name\":" << json::quote(gate_name(g.kind)) << ",\"qubits\":[" << static_cast<unsigned>(g.q[0]);
    if (g.nq == 2) out << ',' << static_cast<unsigned>(g.q[1]);
    out << "],\"params\":[";
    for (int i = 0; i < info(g.kind).nparams; ++i) { if (i) out << ','; out << g.params[i]; }
    out << ']';
    if (g.kind == GateKind::U1Q || g.kind == GateKind::U2Q) {
        out << ",\"matrix\":[";
        for (int i = 0; i < (1 << (2*g.nq)); ++i) {
            if (i) out << ',';
            out << '[' << g.m[i].real() << ',' << g.m[i].imag() << ']';
        }
        out << ']';
    }
    out << '}';
}
void gates_json(std::ostream& out, const std::vector<Gate>& gates) {
    out << '[';
    for (std::size_t i = 0; i < gates.size(); ++i) { if (i) out << ','; gate_json(out,gates[i]); }
    out << ']';
}
Gate parse_gate(const json::Value& v, int n) {
    const auto kind = gate_kind(v.field("name").str());
    const auto& gi = info(kind);
    const bool generic = kind == GateKind::U1Q || kind == GateKind::U2Q;
    if (generic) v.keys({"name","qubits","params","matrix"});
    else v.keys({"name","qubits","params"});
    const auto& qs = v.field("qubits").items();
    const auto& ps = v.field("params").items();
    if (qs.size() != gi.nq || ps.size() != gi.nparams) throw std::invalid_argument("wrong qubit or parameter count");
    std::array<uint8_t,2> q{};
    std::array<double,3> p{};
    for (std::size_t i = 0; i < qs.size(); ++i) q[i] = static_cast<uint8_t>(qs[i].integer(n-1));
    for (std::size_t i = 0; i < ps.size(); ++i) p[i] = ps[i].real();
    if (!generic) return gi.nq == 1 ? make_gate(kind,q[0],p) : make_gate(kind,q[0],q[1]);
    const auto& matrix = v.field("matrix").items();
    if (matrix.size() != static_cast<std::size_t>(1 << (2*gi.nq))) throw std::invalid_argument("wrong matrix dimensions");
    std::array<cplx,16> m{};
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        const auto& z = matrix[i].items();
        if (z.size() != 2) throw std::invalid_argument("complex entry must be [real,imaginary]");
        m[i] = {z[0].real(),z[1].real()};
    }
    if (gi.nq == 1) return make_unitary_1q(q[0],{m[0],m[1],m[2],m[3]});
    return make_unitary_2q(q[0],q[1],m);
}
void metrics_json(std::ostream& out, const PrefixMetrics& m) {
    out << "{\"common_prefix_gates\":" << m.gates << ",\"common_prefix_fraction\":" << m.fraction
        << ",\"a_gates\":" << m.a_gates << ",\"b_gates\":" << m.b_gates << '}';
}
void config_json(std::ostream& out, const GeneratorConfig& c) {
    out << "{\"qubits\":" << c.qubits << ",\"gates\":" << c.gates << ",\"pairs\":" << c.pairs
        << ",\"prefix_group_size\":" << c.prefix_group_size << ",\"prefix_fraction\":" << c.prefix_fraction
        << ",\"equivalent_fraction\":" << c.equivalent_fraction << ",\"seed\":" << json::quote(std::to_string(c.seed))
        << ",\"gate_set\":" << json::quote(c.gate_set) << ",\"angle_mode\":" << json::quote(c.angle_mode)
        << ",\"verify\":" << json::quote(c.verify) << '}';
}
void case_json(std::ostream& out, const CircuitPair& p, const GeneratorConfig& c, const std::string& extension) {
    const auto& m = p.metadata;
    out << "{\"record_type\":\"case\",\"schema_version\":1,\"case_id\":" << json::quote(m.case_id)
        << ",\"expected_relation\":" << json::quote(relation_name(p.expected))
        << ",\"prefix_group_id\":" << m.prefix_group_id << ",\"prefix_identity\":" << json::quote(m.prefix_identity)
        << ",\"seeds\":{\"master\":" << json::quote(std::to_string(m.master_seed))
        << ",\"prefix\":" << json::quote(std::to_string(m.prefix_seed))
        << ",\"suffix\":" << json::quote(std::to_string(m.suffix_seed))
        << ",\"transform\":" << json::quote(std::to_string(m.transform_seed)) << '}'
        << ",\"requested_prefix_fraction\":" << c.prefix_fraction << ",\"planted_prefix_gates\":" << m.planted_prefix_gates
        << ",\"sharing\":{\"raw\":";
    metrics_json(out,m.sharing.raw); out << ",\"normalized\":"; metrics_json(out,m.sharing.normalized);
    out << "},\"global_phase\":{\"a\":" << p.a.global_phase << ",\"b\":" << p.b.global_phase << '}'
        << ",\"transform\":{\"kind\":" << json::quote(transform_name(m.transform.kind))
        << ",\"location\":" << m.transform.location << ",\"phase_adjustment\":" << m.transform.phase_adjustment
        << ",\"proof\":" << json::quote(m.transform.proof) << ",\"removed\":";
    gates_json(out,m.transform.removed); out << ",\"inserted\":"; gates_json(out,m.transform.inserted); out << '}'
        << ",\"negative_mutation\":";
    if (m.negative_qubit < 0) out << "null";
    else out << "{\"kind\":\"append_x\",\"qubit\":" << m.negative_qubit << ",\"location\":" << p.b.size()-1
             << ",\"proof\":\"B = X(q) V, V is equivalent to A; non-scalar X cannot be a global phase.\"}";
    out << ",\"verification\":{\"status\":" << json::quote(m.verification) << ",\"method\":";
    if (m.verification == "passed") out << "\"dense_full_operator_one_global_phase_and_exact_normalization\",\"tolerance\":" << kVerificationTolerance;
    else out << "null,\"tolerance\":null";
    out << "},\"files\":{\"a\":" << json::quote("cases/"+m.case_id+".a"+extension)
        << ",\"b\":" << json::quote("cases/"+m.case_id+".b"+extension) << "}}\n";
}
void close_checked(std::ofstream& out) {
    out.close(); if (!out) throw std::runtime_error("failed to finish output file");
}
} // namespace
void write_ir(std::ostream& out, const Circuit& c) {
    validate_circuit(c); json::numbers(out);
    out << "{\"schema_version\":1,\"n\":" << c.n << ",\"global_phase\":" << c.global_phase << ",\"gates\":";
    gates_json(out,c.gates); out << "}\n";
    if (!out) throw std::runtime_error("failed to write IR");
}
Circuit read_ir(std::istream& in) {
    const auto root = json::Parser(in).parse();
    root.keys({"schema_version","n","global_phase","gates"});
    if (root.field("schema_version").integer(1) != 1) throw std::invalid_argument("unsupported IR schema version");
    Circuit out(root.field("n").integer(kMaxQubits));
    if (out.n < 1) throw std::invalid_argument("IR n must be positive");
    out.global_phase = root.field("global_phase").real();
    for (const auto& v : root.field("gates").items()) out.add(parse_gate(v,out.n));
    validate_circuit(out); return out;
}
void validate_output_options(const OutputOptions& o) {
    if (o.format != "ir" && o.format != "qasm") throw std::invalid_argument("--format must be qasm or ir");
    if (o.destination.empty()) throw std::invalid_argument("--output DIR is required");
    if (fs::exists(o.destination) && (!fs::is_directory(o.destination) || !fs::is_empty(o.destination)))
        throw std::invalid_argument("output destination must be absent or an empty directory");
}
void write_dataset(const OutputOptions& o, const GeneratorConfig& c, const std::vector<CircuitPair>& cases) {
    validate_config(c); validate_output_options(o);
    if (cases.size() != c.pairs) throw std::invalid_argument("case count disagrees with config");
    std::set<std::string> ids;
    // Check circuits and exportability before touching the destination.
    for (const auto& p : cases) {
        const auto& name = p.metadata.case_id;
        if (name.empty() || name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_") != std::string::npos || !ids.insert(name).second)
            throw std::invalid_argument("unsafe or duplicate case ID");
        for (const auto* circuit : {&p.a,&p.b}) {
            validate_circuit(*circuit);
            if (o.format == "qasm") {
                for (const auto& g : circuit->gates) if (g.kind == GateKind::U1Q || g.kind == GateKind::U2Q)
                    throw std::invalid_argument("QASM does not support generic matrices; use --format ir");
            }
        }
    }
    fs::create_directories(o.destination/"cases");
    const std::string ext = o.format == "ir" ? ".ir.json" : ".qasm";
    std::ofstream manifest((o.destination/"manifest.jsonl.partial").string(),std::ios::binary);
    if (!manifest) throw std::runtime_error("cannot create manifest");
    json::numbers(manifest);
    manifest << "{\"record_type\":\"header\",\"schema_version\":1,\"generator_version\":" << json::quote(kGeneratorVersion)
        << ",\"source_commit\":" << json::quote(QCANOPY_SOURCE_COMMIT) << ",\"source_commit_is_build_base\":true"
        << ",\"format\":" << json::quote(o.format) << ",\"dialect\":" << json::quote(o.format == "qasm" ? "OpenQASM 3.0 + stdgates.inc" : "qcanopy-circuit-json-v1")
        << ",\"dense_qubit_cap\":" << kDenseQubitCap << ",\"prefix_tolerance\":" << kPrefixTolerance
        << ",\"equivalent_count_rule\":\"floor(E*B+0.5)\",\"prng\":\"SplitMix64-v1\",\"config\":";
    config_json(manifest,c); manifest << "}\n";
    for (const auto& p : cases) {
        for (int side = 0; side < 2; ++side) {
            const auto path = o.destination/"cases"/(p.metadata.case_id+(side ? ".b" : ".a")+ext);
            if (fs::exists(path)) throw std::runtime_error("refusing to overwrite circuit file");
            std::ofstream out(path.string(),std::ios::binary);
            if (!out) throw std::runtime_error("cannot create circuit file");
            if (o.format == "ir") write_ir(out,side ? p.b : p.a); else write_qasm(out,side ? p.b : p.a);
            close_checked(out);
        }
        case_json(manifest,p,c,ext);
    }
    close_checked(manifest);
    fs::rename(o.destination/"manifest.jsonl.partial",o.destination/"manifest.jsonl");
}
} // namespace qcanopy::testgen
