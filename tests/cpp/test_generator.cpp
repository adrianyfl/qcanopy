#include "qcanopy/generator.hpp"
#include "qcanopy/normalize.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace qcanopy;
using namespace qcanopy::testgen;
namespace {
int failures = 0;
#define CHECK(c) do { if (!(c)) { std::cerr << __LINE__ << ": " #c "\n"; ++failures; } } while (false)
template<class F> void rejects(F f) { try { f(); CHECK(false); } catch (const std::exception&) {} }
bool identical(const Circuit& a, const Circuit& b) {
    if (a.n != b.n || a.size() != b.size() || a.global_phase != b.global_phase) return false;
    for (std::size_t i=0;i<a.size();++i) if (!same_gate(a.gates[i],b.gates[i],0) || a.gates[i].params != b.gates[i].params) return false;
    return true;
}
void transforms() {
    Circuit c(2);
    c.add(make_gate(GateKind::H,1)).add(make_gate(GateKind::X,0)).add(make_gate(GateKind::S,1)).add(make_gate(GateKind::SWAP,0,1));
    c.global_phase = .37;
    for (int k=0;k<7;++k) {
        Circuit b=c;
        const auto t=apply_equivalent_transform(b,static_cast<TransformKind>(k),0,759);
        CHECK(t.kind==static_cast<TransformKind>(k));
        CHECK(equivalent_up_to_phase(c,b));
    }
    Circuit identity(2), concat=c;
    const auto inv=inverse_circuit(c);
    concat.gates.insert(concat.gates.end(),inv.gates.begin(),inv.gates.end()); concat.global_phase+=inv.global_phase;
    CHECK(equivalent_up_to_phase(identity,concat));
    const auto u=reference_unitary(concat), eye=reference_unitary(identity);
    for (std::size_t i=0;i<u.size();++i) CHECK(std::abs(u[i]-eye[i])<1e-12);
    for (int k=0;k<7;++k) {
        auto b=c;
        (void)apply_equivalent_transform(b,static_cast<TransformKind>(k),c.size(),123);
        CHECK(common_prefix(c,b).gates==c.size());
    }
}
void reference_conventions() {
    Circuit a(2), b(2);
    a.add(make_gate(GateKind::CX,0,1)); b.add(make_gate(GateKind::CX,1,0));
    CHECK(!equivalent_up_to_phase(a,b));
    const auto cx=reference_unitary(a);
    CHECK(std::abs(cx[1*4+3]-cplx(1))<1e-12); // |01> -> |11>: q0 control
    CHECK(std::abs(cx[2*4+2]-cplx(1))<1e-12);
    Circuit x(2); x.add(make_gate(GateKind::X,0));
    const auto ux=reference_unitary(x); CHECK(ux[0*4+1]==cplx(1)); CHECK(ux[2*4+3]==cplx(1));
    Circuit i(1), z(1); z.add(make_gate(GateKind::Z,0));
    CHECK(!equivalent_up_to_phase(i,z)); // equal probabilities and column-dependent phases are insufficient
    i.global_phase=.81; Circuit plain(1); CHECK(equivalent_up_to_phase(i,plain));
}
void batches() {
    const auto golden=generate_pair(GeneratorConfig{},0);
    CHECK(golden.metadata.prefix_seed==UINT64_C(6242409861929405762));
    CHECK(golden.metadata.suffix_seed==UINT64_C(8771501796628513655));
    CHECK(golden.metadata.transform_seed==UINT64_C(918932798763059307));
    GeneratorConfig c; c.qubits=3; c.gates=40; c.pairs=21; c.prefix_group_size=4; c.equivalent_fraction=.5;
    const auto batch=generate_batch(c), again=generate_batch(c);
    std::size_t eq=0;
    bool collapsed=false, distinct=false;
    for (std::size_t i=0;i<batch.size();++i) {
        const auto& p=batch[i]; eq+=p.expected==ExpectedRelation::Equivalent;
        CHECK(identical(p.a,again[i].a)); CHECK(identical(p.b,again[i].b));
        CHECK(p.metadata.verification=="passed");
        const auto single=generate_pair(c,i); CHECK(identical(p.a,single.a)); CHECK(identical(p.b,single.b));
        CHECK(p.metadata.sharing.raw.gates>=20);
        CHECK(p.metadata.sharing.normalized.gates==measure_sharing(p.a,p.b).normalized.gates);
        const auto& group=batch[(i/4)*4];
        for (std::size_t j=0;j<20;++j) CHECK(same_gate(p.a.gates[j],group.a.gates[j],0));
        CHECK(p.metadata.prefix_group_id==i/4);
        if (p.expected==ExpectedRelation::Equivalent) {
            const auto na=normalized_copy(p.a), nb=normalized_copy(p.b);
            if (common_prefix(na,nb).gates==na.size() && na.size()==nb.size()) collapsed=true;
            else distinct=true;
        }
    }
    CHECK(eq==11); CHECK(collapsed); CHECK(distinct);
    CHECK(batch.back().metadata.prefix_group_id==5); // partial final group
    c.seed++; CHECK(!identical(batch[0].a,generate_pair(c,0).a));
    for (auto set : {"clifford","clifford-t","parameterized","mixed"}) for (auto mode : {"discrete","continuous","mixed"}) {
        c.gate_set=set; c.angle_mode=mode; c.qubits=1; c.gates=12; c.pairs=7;
        for (const auto& p : generate_batch(c)) { validate_circuit(p.a); CHECK(p.metadata.verification=="passed"); }
    }
}
void boundaries() {
    GeneratorConfig c; c.pairs=7; c.qubits=1;
    for (std::size_t length : {0,1,2}) for (double r : {0.,1.}) for (double e : {0.,1.}) {
        c.gates=length; c.prefix_fraction=r; c.equivalent_fraction=e;
        for (const auto& p : generate_batch(c)) {
            CHECK(p.metadata.planted_prefix_gates==(r==1 ? length:0));
            if (r==1 && !p.a.gates.empty()) CHECK(p.metadata.sharing.raw.fraction==1);
        }
    }
    Circuit empty(1), one(1); one.add(make_gate(GateKind::X,0));
    CHECK(common_prefix(empty,empty).fraction==1); CHECK(common_prefix(empty,one).fraction==0);
    c.qubits=kMaxQubits; c.gates=3; c.pairs=2;
    CHECK(generate_pair(c,0).metadata.verification=="skipped_size_limit");
    rejects([&]{reference_unitary(Circuit(kMaxQubits));});
    c.verify="none"; CHECK(generate_pair(c,0).metadata.verification=="skipped_by_request");
    { auto cap=c; cap.qubits=5; cap.verify="small"; CHECK(generate_pair(cap,0).metadata.verification=="passed"); }
    { auto accidental=c; accidental.qubits=1; accidental.prefix_fraction=0; accidental.equivalent_fraction=1;
      const auto p=generate_pair(accidental,0); CHECK(p.metadata.planted_prefix_gates==0); CHECK(p.metadata.sharing.raw.gates==3); }
    auto bad=c; bad.qubits=0; rejects([&]{validate_config(bad);});
    bad=c; bad.qubits=kMaxQubits+1; rejects([&]{validate_config(bad);});
    bad=c; bad.pairs=0; rejects([&]{validate_config(bad);});
    bad=c; bad.prefix_group_size=0; rejects([&]{validate_config(bad);});
    bad=c; bad.prefix_fraction=std::numeric_limits<double>::quiet_NaN(); rejects([&]{validate_config(bad);});
    bad=c; bad.equivalent_fraction=1.01; rejects([&]{validate_config(bad);});
    bad=c; bad.gates=std::numeric_limits<std::size_t>::max(); rejects([&]{validate_config(bad);});
    Circuit malformed(1); malformed.add(make_gate(GateKind::RX,0,{.2,0,0})); malformed.gates[0].params[0]=.3;
    rejects([&]{validate_circuit(malformed);});
    malformed.gates[0].params[0]=std::numeric_limits<double>::infinity(); rejects([&]{validate_circuit(malformed);});
    malformed.gates[0].kind=GateKind::U1Q; malformed.gates[0].params={}; malformed.gates[0].m={};
    rejects([&]{validate_circuit(malformed);});
    // Explicit normalization tolerance boundaries, kept out of ordinary random angles.
    for (double theta : {.5*kAngleTol,2*kAngleTol}) {
        Circuit near(1); near.add(make_gate(GateKind::P,0,{theta,0,0}));
        CHECK(normalized_copy(near).size()==(theta<kAngleTol ? 0u:1u));
    }
    Circuit cancel(1); cancel.add(make_gate(GateKind::T,0)).add(make_gate(GateKind::Tdg,0));
    CHECK(normalized_copy(cancel).size()==0);
}
}
int main() {
    try { transforms(); reference_conventions(); batches(); boundaries(); }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
    return failures ? 1:0;
}
