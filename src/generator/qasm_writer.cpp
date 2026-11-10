#include "qcanopy/dataset_io.hpp"
#include "json.hpp"
#include <ostream>
#include <stdexcept>

namespace qcanopy::testgen {
void write_qasm(std::ostream& out, const Circuit& c) {
    validate_circuit(c);
    constexpr double pi = 3.14159265358979323846;
    double phase = c.global_phase;
    for (const auto& g : c.gates) {
        if (g.kind == GateKind::U1Q || g.kind == GateKind::U2Q)
            throw std::invalid_argument("OpenQASM 3 export does not support generic U1Q/U2Q matrices; use --format ir");
        if (g.kind == GateKind::SX) phase += pi/4;
        if (g.kind == GateKind::SXdg) phase -= pi/4;
        // OpenQASM 3 U = exp(i*theta/2) * the phase-free IR U.
        if (g.kind == GateKind::U) phase -= g.params[0]/2;
    }
    if (!std::isfinite(phase)) throw std::invalid_argument("exported global phase overflow");
    json::numbers(out);
    out << "OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[" << c.n << "] q;\ngphase(" << phase << ");\n";
    for (const auto& g : c.gates) {
        const auto& gi = info(g.kind);
        if (g.kind == GateKind::SX || g.kind == GateKind::SXdg)
            out << "rx(" << (g.kind == GateKind::SX ? pi/2 : -pi/2) << ')';
        else {
            out << (g.kind == GateKind::U ? "U" : gi.name);
            if (gi.nparams) {
                out << '(';
                for (int i = 0; i < gi.nparams; ++i) { if (i) out << ','; out << g.params[i]; }
                out << ')';
            }
        }
        out << " q[" << static_cast<unsigned>(g.q[0]) << ']';
        if (g.nq == 2) out << ",q[" << static_cast<unsigned>(g.q[1]) << ']';
        out << ";\n";
    }
    if (!out) throw std::runtime_error("failed to write QASM");
}
} // namespace qcanopy::testgen
