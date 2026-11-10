#include "qcanopy/dataset_io.hpp"
#include "../../src/generator/json.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace qcanopy;
using namespace qcanopy::testgen;
namespace {
int failures=0;
#define CHECK(c) do { if (!(c)) { std::cerr << __LINE__ << ": " #c "\n"; ++failures; } } while(false)
template<class F> void rejects(F f) { try { f(); CHECK(false); } catch (const std::exception&) {} }
std::string serialize(const Circuit& c) { std::ostringstream s; write_ir(s,c); return s.str(); }
void roundtrip(const Circuit& c) {
    const auto text=serialize(c); std::istringstream in(text); const auto out=read_ir(in);
    CHECK(out.n==c.n); CHECK(out.global_phase==c.global_phase); CHECK(out.size()==c.size());
    for (std::size_t i=0;i<c.size();++i) {
        CHECK(same_gate(out.gates[i],c.gates[i],0)); CHECK(out.gates[i].params==c.gates[i].params);
        CHECK(out.gates[i].q==c.gates[i].q);
    }
    CHECK(serialize(out)==text);
}
void malformed() {
    const std::string base="{\"schema_version\":1,\"n\":2,\"global_phase\":0,\"gates\":[]}";
    for (auto text : {"", "[]", "{}", "null", "{\"n\":1,\"n\":2}", "{\"schema_version\":1,}", "{\"x\":01}", "{\"x\":1.}", "{\"x\":1e}", "{\"x\":1e999}", "{\"x\":NaN}", "{\"x\":\"\\uD800\"}"})
        rejects([&]{std::istringstream in(text); (void)read_ir(in);});
    auto replace=[&](std::string from,std::string to) {
        auto text=base; text.replace(text.find(from),from.size(),to);
        rejects([&]{std::istringstream in(text); (void)read_ir(in);});
    };
    replace("\"schema_version\":1","\"schema_version\":2");
    replace("\"n\":2","\"n\":0"); replace("\"n\":2","\"n\":33"); replace("\"n\":2","\"n\":1.5");
    replace("\"global_phase\":0","\"global_phase\":1e999");
    for (auto gate : {
        R"({"name":"unknown","qubits":[0],"params":[]})",
        R"({"name":"h","qubits":[2],"params":[]})",
        R"({"name":"rx","qubits":[0],"params":[]})",
        R"({"name":"h","qubits":[0],"params":[1]})",
        R"({"name":"cx","qubits":[0,0],"params":[]})",
        R"({"name":"h","qubits":[0,1],"params":[]})",
        R"({"name":"unitary_1q","qubits":[0],"params":[],"matrix":[[1,0]]})",
        R"({"name":"unitary_1q","qubits":[0],"params":[],"matrix":[[0,0],[0,0],[0,0],[0,0]]})",
        R"({"name":"h","qubits":[0],"params":[],"matrix":[]})"
    }) replace("[]", "["+std::string(gate)+"]");
    rejects([&]{std::istringstream in(base+"garbage"); (void)read_ir(in);});
}
void dataset() {
    const auto root=fs::temp_directory_path()/("qcanopy-io-"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    GeneratorConfig c; c.pairs=9; c.gates=12; c.qubits=2; c.seed=UINT64_MAX;
    const auto cases=generate_batch(c);
    write_dataset({"ir",root/"ir"},c,cases); write_dataset({"qasm",root/"qasm"},c,cases);
    write_dataset({"ir",root/"repeat"},c,generate_batch(c));
    const auto read_text=[](const fs::path& path) { std::ifstream f(path.string(),std::ios::binary); return std::string(std::istreambuf_iterator<char>(f),{}); };
    CHECK(read_text(root/"ir"/"manifest.jsonl")==read_text(root/"repeat"/"manifest.jsonl"));
    for (const auto& pair : cases) for (const char* side : {".a.ir.json",".b.ir.json"})
        CHECK(read_text(root/"ir"/"cases"/(pair.metadata.case_id+side))==read_text(root/"repeat"/"cases"/(pair.metadata.case_id+side)));
    std::ifstream ir((root/"ir"/"manifest.jsonl").string()), qasm((root/"qasm"/"manifest.jsonl").string());
    std::string a,b; std::size_t count=0;
    while (std::getline(ir,a)) {
        CHECK(static_cast<bool>(std::getline(qasm,b)));
        std::istringstream sa(a),sb(b); const auto ja=json::Parser(sa).parse(), jb=json::Parser(sb).parse();
        if (count==0) {
            CHECK(ja.field("format").str()=="ir"); CHECK(jb.field("format").str()=="qasm");
            CHECK(ja.field("config").field("seed").str()=="18446744073709551615");
        } else {
            const auto& p=cases[count-1]; CHECK(ja.field("case_id").str()==p.metadata.case_id);
            CHECK(ja.field("expected_relation").str()==jb.field("expected_relation").str());
            CHECK(ja.field("seeds").field("suffix").str()==jb.field("seeds").field("suffix").str());
            CHECK(ja.field("verification").field("status").str()=="passed");
            std::ifstream file((root/"ir"/ja.field("files").field("a").str()).string());
            CHECK(serialize(read_ir(file))==serialize(p.a));
            CHECK(fs::exists(root/"qasm"/jb.field("files").field("b").str()));
        }
        ++count;
    }
    CHECK(count==10); CHECK(!std::getline(qasm,b));
    ir.close(); qasm.close();
    rejects([&]{write_dataset({"ir",root/"ir"},c,cases);});
    rejects([&]{write_dataset({"ir",root/"ir"/"manifest.jsonl"},c,cases);});
    rejects([&]{write_dataset({"bad",root/"bad"},c,cases);}); CHECK(!fs::exists(root/"bad"));
    auto bad=cases; bad[0].metadata.case_id="../escape";
    rejects([&]{write_dataset({"ir",root/"unsafe"},c,bad);}); CHECK(!fs::exists(root/"unsafe"));
    bad=cases; bad[1].metadata.case_id=bad[0].metadata.case_id;
    rejects([&]{write_dataset({"ir",root/"duplicate"},c,bad);});
    bad=cases; const auto h=make_gate(GateKind::H,0);
    bad[0].a.add(make_unitary_1q(0,{h.m[0],h.m[1],h.m[2],h.m[3]}));
    rejects([&]{write_dataset({"qasm",root/"generic"},c,bad);}); CHECK(!fs::exists(root/"generic"));
    // Only remove the unique temporary directory created by this test.
    fs::remove_all(root);
}
}
int main() {
    try {
        Circuit all(2); all.global_phase=.3141592653589793;
        for (int i=0;i<static_cast<int>(GateKind::U1Q);++i) {
            auto k=static_cast<GateKind>(i);
            std::array<double,3> p{}; for (int j=0;j<info(k).nparams;++j) p[j]=.12345678912345678*(j+1);
            all.add(info(k).nq==1 ? make_gate(k,0,p):make_gate(k,1,0));
        }
        const auto h=make_gate(GateKind::H,0), cx=make_gate(GateKind::CX,0,1);
        all.add(make_unitary_1q(1,{h.m[0],h.m[1],h.m[2],h.m[3]})); all.add(make_unitary_2q(1,0,cx.m));
        roundtrip(all); roundtrip(Circuit(32)); malformed(); dataset();
        Circuit bad(1); bad.add(make_gate(GateKind::X,0)); bad.gates[0].m[0]=1;
        rejects([&]{serialize(bad);});
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
    return failures ? 1:0;
}
