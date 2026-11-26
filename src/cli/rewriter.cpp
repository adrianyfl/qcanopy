#include "rewriter.hpp"

#include "json.hpp"                // src/generator: private strict JSON reader
#include "qcanopy/dataset_io.hpp"  // write_ir

#include <array>
#include <cstdlib>
#include <functional>
#include <filesystem>
#include <fstream>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifndef _WIN32
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
extern char** environ;
#endif

#ifndef QCANOPY_REWRITER_DIR
#define QCANOPY_REWRITER_DIR "rewriting"
#endif
#ifndef QCANOPY_PYTHON_DEFAULT
#define QCANOPY_PYTHON_DEFAULT "python3"
#endif

namespace qcanopy::cli {

    namespace fs = std::filesystem;
    namespace json = testgen::json;

    namespace {
        // IR gate kind for a prover op name; U1Q/U2Q have no name and are never returned.
        GateKind kind_for(const std::string& op) {
            for (int i = 0; i < static_cast<int>(GateKind::Count_); ++i) {
                const auto k = static_cast<GateKind>(i);
                if (!info(k).name.empty() && info(k).name == op) return k;
            }
            throw RewriterError("residual gate '" + op + "' has no IR equivalent");
        }

        // A qdag/1 residual as an IR circuit, in a topological order of its per-wire edges
        // (ties broken by node order, so an already ordered DAG is kept as is).
        Circuit residual_to_circuit(const json::Value& dag, int n) {
            if (dag.field("num_qubits").integer(kMaxQubits) != n)
                throw RewriterError("residual has a different qubit count");
            const auto& nodes = dag.field("nodes").items();
            std::vector<Gate> gates;
            std::vector<std::string> ids;
            for (const auto& node : nodes) {
                const std::string& op = node.field("op").str();
                const GateKind kind = kind_for(op);
                const auto& gi = info(kind);
                const auto& qs = node.field("qubits").items();
                const auto& ps = node.field("params").items();
                if (qs.size() != gi.nq || ps.size() != gi.nparams)
                    throw RewriterError("residual gate '" + op + "' has the wrong arity");
                std::array<uint8_t, 2> q{};
                std::array<double, 3> p{};
                for (std::size_t i = 0; i < qs.size(); ++i) q[i] = static_cast<uint8_t>(qs[i].integer(n - 1));
                for (std::size_t i = 0; i < ps.size(); ++i) {
                    if (ps[i].type != json::Value::Number)
                        throw RewriterError("residual gate '" + op + "' has a symbolic parameter");
                    p[i] = ps[i].real();
                }
                gates.push_back(gi.nq == 1 ? make_gate(kind, q[0], p) : make_gate(kind, q[0], q[1]));
                ids.push_back(node.field("id").str());
            }

            auto index_of = [&](const std::string& id) {
                for (std::size_t i = 0; i < ids.size(); ++i) if (ids[i] == id) return i;
                throw RewriterError("residual edge names unknown gate '" + id + "'");
            };
            std::vector<std::vector<std::size_t>> succ(gates.size());
            std::vector<std::size_t> indegree(gates.size(), 0);
            for (const auto& e : dag.field("edges").items()) {
                const std::size_t s = index_of(e.field("source").str()), t = index_of(e.field("target").str());
                succ[s].push_back(t);
                ++indegree[t];
            }
            std::priority_queue<std::size_t, std::vector<std::size_t>, std::greater<>> ready;
            for (std::size_t i = 0; i < gates.size(); ++i) if (!indegree[i]) ready.push(i);
            Circuit c(n);
            while (!ready.empty()) {
                const std::size_t i = ready.top();
                ready.pop();
                c.add(gates[i]);
                for (const std::size_t t : succ[i]) if (!--indegree[t]) ready.push(t);
            }
            if (c.size() != gates.size()) throw RewriterError("residual DAG has a cycle");
            return c;
        }

        std::string read_file(const fs::path& p) {
            std::ifstream in(p, std::ios::binary);
            std::ostringstream s;
            s << in.rdbuf();
            return s.str();
        }

        // Last non-empty line of the child's stderr, which is where run.py puts "error: ...".
        std::string last_line(const std::string& text) {
            std::size_t end = text.find_last_not_of(" \t\r\n");
            if (end == std::string::npos) return "";
            const std::size_t start = text.find_last_of('\n', end);
            return text.substr(start == std::string::npos ? 0 : start + 1, end - (start == std::string::npos ? 0 : start + 1) + 1);
        }

#ifndef _WIN32
        struct TempDir {
            fs::path path;
            TempDir() {
                std::string tmpl = (fs::temp_directory_path() / "qcanopy-check-XXXXXX").string();
                if (!mkdtemp(tmpl.data())) throw RewriterError(std::string("cannot create temp dir: ") + std::strerror(errno));
                path = tmpl;
            }
            ~TempDir() { std::error_code ec; fs::remove_all(path, ec); }
            TempDir(const TempDir&) = delete;
            TempDir& operator=(const TempDir&) = delete;
        };

        // Runs args (no shell) with stdout discarded and stderr sent to err_file; returns the
        // exit status. Throws RewriterError if the program cannot be started or is killed.
        int spawn_and_wait(const std::vector<std::string>& args, const fs::path& err_file) {
            posix_spawn_file_actions_t actions;
            posix_spawn_file_actions_init(&actions);
            const std::string err_path = err_file.string();
            posix_spawn_file_actions_addopen(&actions, 1, "/dev/null", O_WRONLY, 0);
            posix_spawn_file_actions_addopen(&actions, 2, err_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);

            std::vector<char*> argv;
            for (const auto& a : args) argv.push_back(const_cast<char*>(a.c_str()));
            argv.push_back(nullptr);

            pid_t pid = 0;
            const int rc = posix_spawnp(&pid, argv[0], &actions, nullptr, argv.data(), environ);
            posix_spawn_file_actions_destroy(&actions);
            if (rc != 0) throw RewriterError("cannot run " + args[0] + ": " + std::strerror(rc));

            int status = 0;
            while (waitpid(pid, &status, 0) < 0)
                if (errno != EINTR) throw RewriterError(std::string("waitpid failed: ") + std::strerror(errno));
            if (!WIFEXITED(status)) throw RewriterError("rewriter was killed by signal " + std::to_string(WTERMSIG(status)));
            return WEXITSTATUS(status);
        }
#endif
    } // namespace

    RewriterConfig default_rewriter_config() {
        RewriterConfig cfg;
        const char* py = std::getenv("QCANOPY_PYTHON");
        const char* dir = std::getenv("QCANOPY_REWRITER");
        cfg.python = py && *py ? py : QCANOPY_PYTHON_DEFAULT;
        cfg.dir = dir && *dir ? dir : QCANOPY_REWRITER_DIR;
        return cfg;
    }

    PairProof parse_pair_proof(std::istream& in, int n) {
        try {
            const json::Value root = json::Parser(in).parse();
            if (root.field("schema").str() != "qrewrite-pair/1")
                throw RewriterError("unexpected output schema " + root.field("schema").str());
            const json::Value& proof = root.field("proof");
            const std::string& verdict = proof.field("verdict").str();
            if (verdict != "equivalent_up_to_global_phase" && verdict != "inconclusive")
                throw RewriterError("unexpected verdict " + verdict);

            const json::Value& cores = root.field("layout").field("cores_used");
            if (cores.type != json::Value::Boolean) throw RewriterError("layout.cores_used is not a boolean");

            PairProof out;
            out.proved = verdict == "equivalent_up_to_global_phase";
            out.cores_used = cores.number != 0;
            out.steps = proof.field("proof").items().size();
            out.left = residual_to_circuit(proof.field("residual_left"), n);
            out.right = residual_to_circuit(proof.field("residual_right"), n);
            if (out.proved && (out.left.size() || out.right.size()))
                throw RewriterError("proved verdict with a non-empty residual");
            return out;
        } catch (const RewriterError&) {
            throw;
        } catch (const std::exception& e) {
            throw RewriterError(std::string("unreadable rewriter output: ") + e.what());
        }
    }

    PairProof prove_pair(const Circuit& a, const Circuit& b, const RewriterConfig& cfg) {
        if (a.n != b.n) throw RewriterError("circuits have different qubit counts");
#ifdef _WIN32
        (void)cfg;
        throw RewriterError("running the rewriter is only implemented on POSIX systems");
#else
        const fs::path run_py = fs::path(cfg.dir) / "run.py";
        if (!fs::is_regular_file(run_py)) throw RewriterError("rewriter not found at " + run_py.string());

        TempDir tmp;
        const fs::path left = tmp.path / "a.json", right = tmp.path / "b.json";
        const fs::path out = tmp.path / "pair.json", err = tmp.path / "stderr.txt";
        for (const auto& [path, circuit] : {std::pair{left, &a}, std::pair{right, &b}}) {
            std::ofstream f(path);
            testgen::write_ir(f, *circuit);
            f.close();
            if (!f) throw RewriterError("cannot write " + path.string());
        }

        std::vector<std::string> args{cfg.python, run_py.string(), "pair",
                                      "--left", left.string(), "--right", right.string(),
                                      "--input-format", "qcanopy-ir", "-o", out.string()};
        if (cfg.overlap) args.push_back("--enable-overlap");

        // run.py: 0 proved, 2 inconclusive, 1 input/IO error.
        const int status = spawn_and_wait(args, err);
        if (status != 0 && status != 2) {
            const std::string why = last_line(read_file(err));
            throw RewriterError("rewriter exited with status " + std::to_string(status) + (why.empty() ? "" : ": " + why));
        }
        std::ifstream result(out, std::ios::binary);
        if (!result) throw RewriterError("rewriter produced no output");
        PairProof proof = parse_pair_proof(result, a.n);
        if (proof.proved != (status == 0)) throw RewriterError("rewriter exit status disagrees with its verdict");
        return proof;
#endif
    }

} // namespace qcanopy::cli
