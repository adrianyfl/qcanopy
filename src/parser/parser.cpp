#include "qcanopy/parser.hpp"

#include "qcanopy/lexer.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace qcanopy {

    namespace {
        using qasm::Lexer;
        using qasm::Token;
        using qasm::TokenKind;

        struct GateSpec {
            std::string_view name;
            GateKind kind;
            bool builtin;       // defined by the language itself, not by an included library
        };

        // OpenQASM 3.0: U is built in, the rest come from stdgates.inc.
        constexpr GateSpec kGates3[] = {
            {"U",    GateKind::U,    true},
            {"id",   GateKind::I,    false},
            {"h",    GateKind::H,    false},
            {"x",    GateKind::X,    false},
            {"y",    GateKind::Y,    false},
            {"z",    GateKind::Z,    false},
            {"s",    GateKind::S,    false},
            {"sdg",  GateKind::Sdg,  false},
            {"t",    GateKind::T,    false},
            {"tdg",  GateKind::Tdg,  false},
            {"rx",   GateKind::RX,   false},
            {"ry",   GateKind::RY,   false},
            {"rz",   GateKind::RZ,   false},
            {"p",    GateKind::P,    false},
            {"cx",   GateKind::CX,   false},
            {"cz",   GateKind::CZ,   false},
            {"swap", GateKind::SWAP, false},
        };

        // OpenQASM 2.0: U and CX are built in; the rest are the qelib1.inc gates the IR has.
        constexpr GateSpec kGates2[] = {
            {"U",    GateKind::U,    true},
            {"CX",   GateKind::CX,   true},
            {"id",   GateKind::I,    false},
            {"x",    GateKind::X,    false},
            {"y",    GateKind::Y,    false},
            {"z",    GateKind::Z,    false},
            {"h",    GateKind::H,    false},
            {"s",    GateKind::S,    false},
            {"sdg",  GateKind::Sdg,  false},
            {"t",    GateKind::T,    false},
            {"tdg",  GateKind::Tdg,  false},
            {"rx",   GateKind::RX,   false},
            {"rz",   GateKind::RZ,   false},
            {"cx",   GateKind::CX,   false},
            {"cz",   GateKind::CZ,   false},
        };

        // qelib1.inc from the OpenQASM 2 specification (github.com/openqasm/openqasm, branch
        // OpenQASM2.x; Apache-2.0), plus cu3's "u1((lambda+phi)/2) c;" (parser.hpp, point 3).
        constexpr std::string_view kQelib1 = R"qelib(// Quantum Experience (QE) Standard Header
// file: qelib1.inc

// --- QE Hardware primitives ---

// 3-parameter 2-pulse single qubit gate
gate u3(theta,phi,lambda) q { U(theta,phi,lambda) q; }
// 2-parameter 1-pulse single qubit gate
gate u2(phi,lambda) q { U(pi/2,phi,lambda) q; }
// 1-parameter 0-pulse single qubit gate
gate u1(lambda) q { U(0,0,lambda) q; }
// controlled-NOT
gate cx c,t { CX c,t; }
// idle gate (identity)
gate id a { U(0,0,0) a; }

// --- QE Standard Gates ---

// Pauli gate: bit-flip
gate x a { u3(pi,0,pi) a; }
// Pauli gate: bit and phase flip
gate y a { u3(pi,pi/2,pi/2) a; }
// Pauli gate: phase flip
gate z a { u1(pi) a; }
// Clifford gate: Hadamard
gate h a { u2(0,pi) a; }
// Clifford gate: sqrt(Z) phase gate
gate s a { u1(pi/2) a; }
// Clifford gate: conjugate of sqrt(Z)
gate sdg a { u1(-pi/2) a; }
// C3 gate: sqrt(S) phase gate
gate t a { u1(pi/4) a; }
// C3 gate: conjugate of sqrt(S)
gate tdg a { u1(-pi/4) a; }

// --- Standard rotations ---
// Rotation around X-axis
gate rx(theta) a { u3(theta,-pi/2,pi/2) a; }
// rotation around Y-axis
gate ry(theta) a { u3(theta,0,0) a; }
// rotation around Z axis
gate rz(phi) a { u1(phi) a; }

// --- QE Standard User-Defined Gates  ---

// controlled-Phase
gate cz a,b { h b; cx a,b; h b; }
// controlled-Y
gate cy a,b { sdg b; cx a,b; s b; }
// controlled-H
gate ch a,b {
h b; sdg b;
cx a,b;
h b; t b;
cx a,b;
t b; h b; s b; x b; s a;
}
// C3 gate: Toffoli
gate ccx a,b,c
{
  h c;
  cx b,c; tdg c;
  cx a,c; t c;
  cx b,c; tdg c;
  cx a,c; t b; t c; h c;
  cx a,b; t a; tdg b;
  cx a,b;
}
// controlled rz rotation
gate crz(lambda) a,b
{
  u1(lambda/2) b;
  cx a,b;
  u1(-lambda/2) b;
  cx a,b;
}
// controlled phase rotation
gate cu1(lambda) a,b
{
  u1(lambda/2) a;
  cx a,b;
  u1(-lambda/2) b;
  cx a,b;
  u1(lambda/2) b;
}
// controlled-U
gate cu3(theta,phi,lambda) c, t
{
  // implements controlled-U(theta,phi,lambda) with  target t and control c
  u1((lambda+phi)/2) c;
  u1((lambda-phi)/2) t;
  cx c,t;
  u3(-theta/2,0,-(phi+lambda)/2) t;
  cx c,t;
  u3(theta/2,phi,0) t;
}
)qelib";

        template <std::size_t N>
        const GateSpec* find_in(const GateSpec (&table)[N], std::string_view name) {
            for (const GateSpec& g : table)
                if (g.name == name) return &g;
            return nullptr;
        }

        // A qubit register: qubits offset .. offset + size - 1 of the circuit.
        struct Register {
            std::string_view name;
            int offset;
            int size;
        };

        // Expression node. A tree, so a gate body can be evaluated again at each call.
        struct Expr {
            enum class Op : uint8_t { Number, Param, Neg, Add, Sub, Mul, Div };
            Op op;
            double value;       // Number
            int a, b;           // Param: index a; Neg: operand a; Add .. Div: operands a, b
        };

        double evaluate(const std::vector<Expr>& pool, int i, const std::vector<double>& params) {
            const Expr& e = pool[i];
            switch (e.op) {
                case Expr::Op::Number: return e.value;
                case Expr::Op::Param:  return params[e.a];
                case Expr::Op::Neg:    return -evaluate(pool, e.a, params);
                case Expr::Op::Add:    return evaluate(pool, e.a, params) + evaluate(pool, e.b, params);
                case Expr::Op::Sub:    return evaluate(pool, e.a, params) - evaluate(pool, e.b, params);
                case Expr::Op::Mul:    return evaluate(pool, e.a, params) * evaluate(pool, e.b, params);
                case Expr::Op::Div:    return evaluate(pool, e.a, params) / evaluate(pool, e.b, params);
            }
            return 0;
        }

        // A resolved gate name: a table entry, or an index into defs_.
        struct Target {
            const GateSpec* spec = nullptr;
            int def = -1;

            bool found() const { return spec || def >= 0; }
        };

        // One statement of a gate body.
        struct Call {
            Target target;
            std::vector<int> params;    // expression roots in the definition's pool
            std::vector<int> args;      // indices into the definition's qubit arguments
        };

        struct GateDef {
            std::string_view name;
            int nparams = 0;
            int nqubits = 0;
            std::vector<Expr> pool;
            std::vector<Call> body;
        };

        constexpr double kPi = 3.14159265358979323846;

        std::string quoted(std::string_view s) { return "'" + std::string(s) + "'"; }

        // How a token is named in an error message.
        std::string describe(const Token& t) {
            if (t.kind == TokenKind::End) return "end of input";
            if (t.kind == TokenKind::String) return "string \"" + std::string(t.text) + "\"";
            return quoted(t.text);
        }

        std::string count(int n, const char* noun) {
            return std::to_string(n) + " " + noun + (n == 1 ? "" : "s");
        }

        int index_of(const std::vector<std::string_view>& names, std::string_view name) {
            for (std::size_t i = 0; i < names.size(); ++i)
                if (names[i] == name) return static_cast<int>(i);
            return -1;
        }

        // Recursive descent with one token of lookahead (tok_).
        class Parser {
        public:
            explicit Parser(std::string_view source) : lexer_(source), tok_(lexer_.next()) {}

            Circuit parse() {
                version();
                while (tok_.kind != TokenKind::End) statement();
                if (registers_.empty()) fail(tok_, "no qubit register declared");
                return std::move(circuit_);
            }

        private:
            [[noreturn]] void fail(const Token& at, const std::string& message) const {
                throw ParseError(at.line, at.column, message);
            }

            Token advance() {
                const Token t = tok_;
                tok_ = lexer_.next();
                return t;
            }

            bool accept(TokenKind kind) {
                if (tok_.kind != kind) return false;
                advance();
                return true;
            }

            Token expect(TokenKind kind, const std::string& what) {
                if (tok_.kind != kind)
                    fail(tok_, "expected " + what + " but found " + describe(tok_));
                return advance();
            }

            void version() {
                if (tok_.kind != TokenKind::Identifier || tok_.text != "OPENQASM")
                    fail(tok_, "expected a version statement such as 'OPENQASM 3.0;' but found " +
                               describe(tok_));
                advance();
                if (tok_.kind != TokenKind::Integer && tok_.kind != TokenKind::Float)
                    fail(tok_, "expected a version number but found " + describe(tok_));
                if (tok_.text == "2.0")
                    v2_ = true;
                else if (tok_.text != "3.0" && tok_.text != "3")
                    fail(tok_, "unsupported version " + quoted(tok_.text) + "; only 2.0 and 3.0 are supported");
                advance();
                expect(TokenKind::Semicolon, "';'");
            }

            void statement() {
                if (tok_.kind != TokenKind::Identifier)
                    fail(tok_, "expected a statement but found " + describe(tok_));
                if (tok_.text == "OPENQASM")
                    fail(tok_, "the version statement must come first and appear only once");
                if (tok_.text == "include") return include();
                if (v2_) {
                    if (tok_.text == "qreg") return declaration();
                    if (tok_.text == "creg") return classical_register();
                    if (tok_.text == "gate") return definition();
                } else {
                    if (tok_.text == "qubit") return declaration();
                    if (tok_.text == "gphase") return gphase();
                }
                gate();
            }

            // The one gate library each version may include.
            std::string_view library() const { return v2_ ? "qelib1.inc" : "stdgates.inc"; }

            void include() {
                const Token keyword = advance();
                const Token file = expect(TokenKind::String, "a file name in double quotes");
                if (file.text != library())
                    fail(file, "cannot include \"" + std::string(file.text) + "\"; only \"" +
                               std::string(library()) + "\" is supported");
                expect(TokenKind::Semicolon, "';'");
                if (included_) return;

                if (v2_) {
                    for (const GateSpec& g : kGates2)
                        if (find_definition(g.name) >= 0)
                            fail(keyword, "gate " + quoted(g.name) + " is already defined, so \"qelib1.inc\" "
                                          "cannot be included");
                }
                included_ = true;
                // Earlier definitions may clash with the library's: read it now to report that here.
                if (v2_ && !defs_.empty()) read_library(&keyword);
            }

            // The library's definitions are read on first use; most files never need them.
            void need_library() {
                if (v2_ && included_ && !library_read_) read_library(nullptr);
            }

            // Reads qelib1.inc's definitions, then resumes the file. A clash is reported at include.
            void read_library(const Token* include) {
                library_read_ = true;
                const Lexer saved_lexer = lexer_;
                const Token saved_tok = tok_;
                lexer_ = Lexer(kQelib1);
                tok_ = lexer_.next();
                in_library_ = true;
                include_at_ = include;
                while (tok_.kind != TokenKind::End) definition();
                in_library_ = false;
                include_at_ = nullptr;
                lexer_ = saved_lexer;
                tok_ = saved_tok;
            }

            // "qubit[n] name;" in 3.0, "qreg name[n];" in 2.0.
            void declaration() {
                const Token keyword = advance();
                if (!v2_ && !registers_.empty()) fail(keyword, "only one qubit register is supported");

                Token name, size;
                if (v2_) name = expect(TokenKind::Identifier, "a register name");
                expect(TokenKind::LBracket, "'['");
                size = expect(TokenKind::Integer, "a register size");
                expect(TokenKind::RBracket, "']'");
                if (!v2_) name = expect(TokenKind::Identifier, "a register name");
                expect(TokenKind::Semicolon, "';'");

                const unsigned long long n = integer(size);
                if (n < 1) fail(size, "a register needs at least 1 qubit");
                if (n > static_cast<unsigned long long>(kMaxQubits - circuit_.n))
                    fail(size, "circuit would have more than " + std::to_string(kMaxQubits) + " qubits");
                if (find_register(name.text))
                    fail(name, "register " + quoted(name.text) + " is already declared");

                registers_.push_back(Register{name.text, circuit_.n, static_cast<int>(n)});
                circuit_.n += static_cast<int>(n);
            }

            // "creg name[n];": checked and ignored.
            void classical_register() {
                advance();
                expect(TokenKind::Identifier, "a register name");
                expect(TokenKind::LBracket, "'['");
                expect(TokenKind::Integer, "a register size");
                expect(TokenKind::RBracket, "']'");
                expect(TokenKind::Semicolon, "';'");
            }

            const Register* find_register(std::string_view name) const {
                for (const Register& r : registers_)
                    if (r.name == name) return &r;
                return nullptr;
            }

            void gphase() {
                advance();
                expect(TokenKind::LParen, "'('");
                circuit_.global_phase += argument();
                expect(TokenKind::RParen, "')'");
                expect(TokenKind::Semicolon, "';'");
            }

            // ---- Gate names ----

            const GateSpec* find_spec(std::string_view name) const {
                return v2_ ? find_in(kGates2, name) : find_in(kGates3, name);
            }

            int find_definition(std::string_view name) const {
                for (std::size_t i = 0; i < defs_.size(); ++i)
                    if (defs_[i].name == name) return static_cast<int>(i);
                return -1;
            }

            Target resolve(std::string_view name) const {
                const GateSpec* spec = find_spec(name);
                if (spec && (spec->builtin || included_)) return Target{spec, -1};
                return Target{nullptr, find_definition(name)};
            }

            // The gate a statement names; fails if there is none.
            Target callee(const Token& name) {
                Target t = resolve(name.text);
                if (!t.found()) {
                    need_library();
                    t = resolve(name.text);
                }
                if (t.found()) return t;
                if (find_spec(name.text))
                    fail(name, "gate " + quoted(name.text) + " is not defined; include \"" +
                               std::string(library()) + "\" first");
                fail(name, "unknown gate or statement " + quoted(name.text));
            }

            int nparams_of(Target t) const { return t.spec ? info(t.spec->kind).nparams : defs_[t.def].nparams; }
            int nqubits_of(Target t) const { return t.spec ? info(t.spec->kind).nq : defs_[t.def].nqubits; }

            void check_params(const Token& name, Target t, std::size_t n) const {
                if (static_cast<int>(n) != nparams_of(t))
                    fail(name, "gate " + quoted(name.text) + " takes " + count(nparams_of(t), "parameter") +
                               " but got " + std::to_string(n));
            }

            // ids: the operands; at: where each is written.
            void check_operands(const Token& name, Target t, const std::vector<int>& ids,
                                const std::vector<Token>& at) const {
                if (static_cast<int>(ids.size()) != nqubits_of(t))
                    fail(name, "gate " + quoted(name.text) + " acts on " + count(nqubits_of(t), "qubit") +
                               " but got " + std::to_string(ids.size()));
                for (std::size_t i = 1; i < ids.size(); ++i)
                    for (std::size_t j = 0; j < i; ++j)
                        if (ids[i] == ids[j])
                            fail(at[i], "gate " + quoted(name.text) + " needs distinct qubits");
            }

            // ---- Gate statements ----

            void gate() {
                const Token name = advance();
                const Target target = callee(name);

                // Member buffers, so a gate statement allocates nothing.
                std::vector<double>& params = call_params_;
                std::vector<int>& qubits = call_qubits_;
                std::vector<Token>& at = call_at_;
                params.clear();
                qubits.clear();
                at.clear();

                if (accept(TokenKind::LParen)) {
                    if (tok_.kind != TokenKind::RParen) {
                        do params.push_back(argument());
                        while (accept(TokenKind::Comma));
                    }
                    expect(TokenKind::RParen, "')'");
                }
                check_params(name, target, params.size());

                do {
                    at.push_back(tok_);
                    qubits.push_back(qubit());
                } while (accept(TokenKind::Comma));
                expect(TokenKind::Semicolon, "';'");
                check_operands(name, target, qubits, at);

                apply(name, target, params, qubits);
            }

            // Adds the gates of one call, expanding defined gates.
            void apply(const Token& at, Target t, const std::vector<double>& params, const std::vector<int>& qubits) {
                if (t.spec) {
                    const GateKind kind = t.spec->kind;
                    if (info(kind).nq == 1) {
                        std::array<double, 3> p{0, 0, 0};
                        for (std::size_t i = 0; i < params.size(); ++i) p[i] = params[i];
                        circuit_.add(make_gate(kind, static_cast<uint8_t>(qubits[0]), p));
                    } else {
                        circuit_.add(make_gate(kind, static_cast<uint8_t>(qubits[0]), static_cast<uint8_t>(qubits[1])));
                    }
                    // parser.hpp, point 3.
                    if (!v2_ && kind == GateKind::U) circuit_.global_phase += params[0] / 2;
                    return;
                }

                const GateDef& def = defs_[t.def];
                std::vector<double> values;
                std::vector<int> operands;
                for (const Call& call : def.body) {
                    values.clear();
                    for (int root : call.params) {
                        const double v = evaluate(def.pool, root, params);
                        if (!std::isfinite(v))
                            fail(at, "gate " + quoted(at.text) + ": a parameter of " + quoted(def.name) +
                                     " does not evaluate to a finite number");
                        values.push_back(v);
                    }
                    operands.clear();
                    for (int a : call.args) operands.push_back(qubits[a]);
                    apply(at, call.target, values, operands);
                }
            }

            int qubit() {
                const Token name = expect(TokenKind::Identifier, "a qubit");
                const Register* reg = find_register(name.text);
                if (!reg) fail(name, "undefined register " + quoted(name.text));
                expect(TokenKind::LBracket, "'['");
                const Token index = expect(TokenKind::Integer, "a qubit index");
                const unsigned long long i = integer(index);
                if (i >= static_cast<unsigned long long>(reg->size))
                    fail(index, "qubit index " + std::string(index.text) + " is out of range for register " +
                                quoted(reg->name) + " of size " + std::to_string(reg->size));
                expect(TokenKind::RBracket, "']'");
                return reg->offset + static_cast<int>(i);
            }

            // ---- Gate definitions (OpenQASM 2.0) ----

            // "gate name(p, ...) a, ... { body }"
            void definition() {
                advance();
                if (!in_library_) need_library();       // to see a clash with a library gate
                const Token name = expect(TokenKind::Identifier, "a gate name");

                std::vector<std::string_view> params, args;
                if (accept(TokenKind::LParen)) {
                    if (tok_.kind != TokenKind::RParen) {
                        do params.push_back(new_name(params, "parameter"));
                        while (accept(TokenKind::Comma));
                    }
                    expect(TokenKind::RParen, "')'");
                }
                do args.push_back(new_name(args, "qubit argument"));
                while (accept(TokenKind::Comma));
                expect(TokenKind::LBrace, "'{'");

                GateDef def;
                def.name = name.text;
                def.nparams = static_cast<int>(params.size());
                def.nqubits = static_cast<int>(args.size());
                pool_ = &def.pool;
                params_ = &params;
                while (!accept(TokenKind::RBrace)) def.body.push_back(body_call(args));
                pool_ = &scratch_;
                params_ = nullptr;

                if (in_library_) {
                    // Table gates keep their IR mapping.
                    if (find_spec(name.text)) return;
                    if (find_definition(name.text) >= 0)
                        fail(include_at_ ? *include_at_ : name,
                             "gate " + quoted(name.text) + " is already defined, so \"qelib1.inc\" "
                             "cannot be included");
                } else if (resolve(name.text).found()) {
                    fail(name, "gate " + quoted(name.text) + " is already defined");
                }
                defs_.push_back(std::move(def));
            }

            std::string_view new_name(const std::vector<std::string_view>& seen, const char* what) {
                const Token t = expect(TokenKind::Identifier, std::string("a ") + what + " name");
                if (index_of(seen, t.text) >= 0) fail(t, std::string("duplicate ") + what + " " + quoted(t.text));
                return t.text;
            }

            // One body statement; its operands are argument names.
            Call body_call(const std::vector<std::string_view>& args) {
                const Token name = expect(TokenKind::Identifier, "a gate or '}'");
                Call call;
                call.target = callee(name);

                if (accept(TokenKind::LParen)) {
                    if (tok_.kind != TokenKind::RParen) {
                        do call.params.push_back(sum());
                        while (accept(TokenKind::Comma));
                    }
                    expect(TokenKind::RParen, "')'");
                }
                check_params(name, call.target, call.params.size());

                std::vector<Token> at;
                do {
                    const Token a = expect(TokenKind::Identifier, "a qubit argument");
                    const int i = index_of(args, a.text);
                    if (i < 0) fail(a, "unknown qubit argument " + quoted(a.text));
                    at.push_back(a);
                    call.args.push_back(i);
                } while (accept(TokenKind::Comma));
                expect(TokenKind::Semicolon, "';'");
                check_operands(name, call.target, call.args, at);
                return call;
            }

            // ---- Expressions (parser.hpp, point 4) ----

            // An argument outside a gate body, evaluated now.
            double argument() {
                const Token start = tok_;
                scratch_.clear();
                const double v = evaluate(scratch_, sum(), {});
                if (!std::isfinite(v))
                    fail(start, "expression does not evaluate to a finite number");
                return v;
            }

            // Each returns the root of its expression in *pool_.
            int sum() {
                int lhs = product();
                while (tok_.kind == TokenKind::Plus || tok_.kind == TokenKind::Minus) {
                    const Expr::Op op = advance().kind == TokenKind::Plus ? Expr::Op::Add : Expr::Op::Sub;
                    const int rhs = product();
                    lhs = node(Expr{op, 0, lhs, rhs});
                }
                return lhs;
            }

            int product() {
                int lhs = unary();
                while (tok_.kind == TokenKind::Star || tok_.kind == TokenKind::Slash) {
                    const Expr::Op op = advance().kind == TokenKind::Star ? Expr::Op::Mul : Expr::Op::Div;
                    const int rhs = unary();
                    lhs = node(Expr{op, 0, lhs, rhs});
                }
                return lhs;
            }

            int unary() {
                if (accept(TokenKind::Minus)) return node(Expr{Expr::Op::Neg, 0, unary(), 0});
                return primary();
            }

            int primary() {
                if (accept(TokenKind::LParen)) {
                    const int v = sum();
                    expect(TokenKind::RParen, "')'");
                    return v;
                }
                if (tok_.kind == TokenKind::Integer || tok_.kind == TokenKind::Float)
                    return node(Expr{Expr::Op::Number, real(advance()), 0, 0});
                if (tok_.kind == TokenKind::Identifier) {
                    const Token name = advance();
                    if (params_) {
                        const int i = index_of(*params_, name.text);
                        if (i >= 0) return node(Expr{Expr::Op::Param, 0, i, 0});
                    }
                    if (name.text == "pi") return node(Expr{Expr::Op::Number, kPi, 0, 0});
                    fail(name, std::string(params_ ? "unknown parameter " : "unknown constant ") + quoted(name.text));
                }
                fail(tok_, "expected a number but found " + describe(tok_));
            }

            int node(const Expr& e) {
                pool_->push_back(e);
                return static_cast<int>(pool_->size() - 1);
            }

            double real(const Token& t) const {
                // The classic locale keeps '.' as the decimal point whatever the global locale is.
                std::istringstream in{std::string(t.text)};
                in.imbue(std::locale::classic());
                double v = 0;
                in >> v;
                if (in.fail() || !std::isfinite(v))
                    fail(t, "number " + quoted(t.text) + " is out of range");
                return v;
            }

            unsigned long long integer(const Token& t) {
                unsigned long long v = 0;
                const auto result = std::from_chars(t.text.data(), t.text.data() + t.text.size(), v);
                if (result.ec != std::errc())
                    fail(t, "integer " + quoted(t.text) + " is too large");
                return v;
            }

            Lexer lexer_;
            Token tok_;                         // next unconsumed token
            Circuit circuit_;
            std::vector<Register> registers_;   // in declaration order; circuit_.n is their total
            std::vector<GateDef> defs_;         // from qelib1.inc and the file
            bool v2_ = false;                   // OpenQASM 2.0, not 3.0
            bool included_ = false;             // include of library() seen
            bool library_read_ = false;         // qelib1.inc's definitions are in defs_
            bool in_library_ = false;           // reading qelib1.inc, not the file
            const Token* include_at_ = nullptr; // while in_library_: the include statement, if known

            // Buffers for gate().
            std::vector<double> call_params_;
            std::vector<int> call_qubits_;
            std::vector<Token> call_at_;

            // Expressions go to scratch_, or to a definition's pool (parameter names in params_).
            std::vector<Expr> scratch_;
            std::vector<Expr>* pool_ = &scratch_;
            const std::vector<std::string_view>* params_ = nullptr;
        };
    }

    Circuit parse_qasm(std::string_view source) {
        return Parser(source).parse();
    }

    Circuit parse_qasm_file(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("cannot open " + path);
        std::ostringstream text;
        text << in.rdbuf();
        if (in.bad()) throw std::runtime_error("cannot read " + path);

        try {
            return parse_qasm(text.str());
        } catch (const ParseError& e) {
            throw ParseError(path, e);
        }
    }

    CircuitDAG parse_qasm_dag(std::string_view source) {
        return CircuitDAG(parse_qasm(source));
    }

    CircuitDAG parse_qasm_file_dag(const std::string& path) {
        return CircuitDAG(parse_qasm_file(path));
    }

}
