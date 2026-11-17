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

        // OpenQASM 2.0: gates from qelib1.inc (parser.hpp, point 3).
        constexpr GateSpec kGates2[] = {
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

        // Recursive descent with one token of lookahead: tok_ is the next unconsumed token,
        // and each method consumes exactly the construct it is named after.
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

            Token expect(TokenKind kind, const char* what) {
                if (tok_.kind != kind)
                    fail(tok_, std::string("expected ") + what + " but found " + describe(tok_));
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
                } else {
                    if (tok_.text == "qubit") return declaration();
                    if (tok_.text == "gphase") return gphase();
                }
                gate();
            }

            // The one gate library each version may include.
            std::string_view library() const { return v2_ ? "qelib1.inc" : "stdgates.inc"; }

            void include() {
                advance();
                const Token file = expect(TokenKind::String, "a file name in double quotes");
                if (file.text != library())
                    fail(file, "cannot include \"" + std::string(file.text) + "\"; only \"" +
                               std::string(library()) + "\" is supported");
                expect(TokenKind::Semicolon, "';'");
                included_ = true;
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

            // "creg name[n];": classical bits are not used, so the statement is only checked.
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

            void gate() {
                const Token name = advance();
                const GateSpec* spec = v2_ ? find_in(kGates2, name.text) : find_in(kGates3, name.text);
                if (!spec) fail(name, "unknown gate or statement " + quoted(name.text));
                if (!spec->builtin && !included_)
                    fail(name, "gate " + quoted(name.text) + " is not defined; include \"" +
                               std::string(library()) + "\" first");
                const GateInfo& gi = info(spec->kind);

                std::array<double, 3> params{0, 0, 0};
                int nparams = 0;
                if (accept(TokenKind::LParen)) {
                    if (tok_.kind != TokenKind::RParen) {
                        do {
                            const double v = argument();
                            if (nparams < 3) params[nparams] = v;
                            ++nparams;
                        } while (accept(TokenKind::Comma));
                    }
                    expect(TokenKind::RParen, "')'");
                }
                if (nparams != gi.nparams)
                    fail(name, "gate " + quoted(name.text) + " takes " + count(gi.nparams, "parameter") +
                               " but got " + std::to_string(nparams));

                std::array<uint8_t, 2> q{0, 0};
                int nq = 0;
                Token second = tok_;
                do {
                    if (nq == 1) second = tok_;
                    const uint8_t index = qubit();
                    if (nq < 2) q[nq] = index;
                    ++nq;
                } while (accept(TokenKind::Comma));
                expect(TokenKind::Semicolon, "';'");
                if (nq != gi.nq)
                    fail(name, "gate " + quoted(name.text) + " acts on " + count(gi.nq, "qubit") +
                               " but got " + std::to_string(nq));
                if (nq == 2 && q[0] == q[1])
                    fail(second, "gate " + quoted(name.text) + " needs two distinct qubits");

                if (nq == 1) circuit_.add(make_gate(spec->kind, q[0], params));
                else         circuit_.add(make_gate(spec->kind, q[0], q[1]));

                // parser.hpp, point 3.
                if (spec->kind == GateKind::U) circuit_.global_phase += params[0] / 2;
            }

            uint8_t qubit() {
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
                return static_cast<uint8_t>(reg->offset + static_cast<int>(i));
            }

            // A gate or gphase argument: a constant expression (parser.hpp, point 4).
            double argument() {
                const Token start = tok_;
                const double v = sum();
                if (!std::isfinite(v))
                    fail(start, "expression does not evaluate to a finite number");
                return v;
            }

            double sum() {
                double v = product();
                while (tok_.kind == TokenKind::Plus || tok_.kind == TokenKind::Minus) {
                    const bool plus = advance().kind == TokenKind::Plus;
                    const double rhs = product();
                    v = plus ? v + rhs : v - rhs;
                }
                return v;
            }

            double product() {
                double v = unary();
                while (tok_.kind == TokenKind::Star || tok_.kind == TokenKind::Slash) {
                    const bool times = advance().kind == TokenKind::Star;
                    const double rhs = unary();
                    v = times ? v * rhs : v / rhs;
                }
                return v;
            }

            double unary() {
                if (accept(TokenKind::Minus)) return -unary();
                return primary();
            }

            double primary() {
                if (accept(TokenKind::LParen)) {
                    const double v = sum();
                    expect(TokenKind::RParen, "')'");
                    return v;
                }
                if (tok_.kind == TokenKind::Integer || tok_.kind == TokenKind::Float)
                    return real(advance());
                if (tok_.kind == TokenKind::Identifier) {
                    const Token name = advance();
                    if (name.text == "pi") return kPi;
                    fail(name, "unknown constant " + quoted(name.text));
                }
                fail(tok_, "expected a number but found " + describe(tok_));
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
            Token tok_;                     // next unconsumed token
            Circuit circuit_;
            std::vector<Register> registers_;   // in declaration order; circuit_.n is their total
            bool v2_ = false;                   // the source is OpenQASM 2.0, not 3.0
            bool included_ = false;             // the include of library() has been seen
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
