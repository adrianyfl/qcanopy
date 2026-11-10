#include "qcanopy/dataset_io.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace qcanopy;
using namespace qcanopy::testgen;
int main() {
    Circuit c(2); c.global_phase=.3;
    for (int i=0;i<static_cast<int>(GateKind::U1Q);++i) {
        auto k=static_cast<GateKind>(i); std::array<double,3> p{};
        for (int j=0;j<info(k).nparams;++j) p[j]=.25*(j+1);
        c.add(info(k).nq==1 ? make_gate(k,0,p):make_gate(k,1,0));
    }
    std::ostringstream out; write_qasm(out,c); const auto s=out.str();
    for (const char* expected : {"OPENQASM 3.0;", "include \"stdgates.inc\";", "qubit[2] q;", "id q[0];", "U(0.25,0.5,0.75) q[0];", "cx q[1],q[0];", "rx(1.5707963267948966) q[0];", "rx(-1.5707963267948966) q[0];"})
        if (s.find(expected)==std::string::npos) { std::cerr << "missing: " << expected << '\n'; return 1; }
    if (s.find("sxdg")!=std::string::npos || s.find("\nu(")!=std::string::npos) return 1;
    Circuit u_gate(1); u_gate.add(make_gate(GateKind::U,0,{.5,.2,.3}));
    std::ostringstream u_text; write_qasm(u_text,u_gate);
    if (u_text.str().find("gphase(-0.25);")==std::string::npos) return 1;
    constexpr double pi=3.14159265358979323846;
    for (auto kind : {GateKind::SX,GateKind::SXdg}) {
        Circuit raw(1), lowered(1); raw.add(make_gate(kind,0));
        lowered.add(make_gate(GateKind::RX,0,{kind==GateKind::SX ? pi/2:-pi/2,0,0}));
        lowered.global_phase=kind==GateKind::SX ? pi/4:-pi/4;
        const auto a=reference_unitary(raw),b=reference_unitary(lowered);
        for (std::size_t i=0;i<a.size();++i) if (std::abs(a[i]-b[i])>1e-14) return 1;
        std::ostringstream q; write_qasm(q,raw);
        std::ostringstream phase; phase.precision(17); phase << "gphase(" << lowered.global_phase << ");";
        if (q.str().find(phase.str())==std::string::npos) return 1;
    }
    const auto h=make_gate(GateKind::H,0);
    Circuit generic(1); generic.add(make_unitary_1q(0,{h.m[0],h.m[1],h.m[2],h.m[3]}));
    try { std::ostringstream unused; write_qasm(unused,generic); return 1; } catch (const std::invalid_argument&) {}
    return 0;
}
