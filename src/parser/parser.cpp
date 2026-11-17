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

namespace qcanopy {

    namespace {
        using qasm::Lexer;
        using qasm::Token;
        using qasm::TokenKind;

        struct GateSpec {
            std::string_view name;
            GateKind kind;
            bool builtin;       // defined by the language itself, not by stdgates.inc
        };

        constexpr GateSpec kGates[] = {
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

        const GateSpec* find_gate(std::string_view name) {
            for (const GateSpec& g : kGates)
                if (g.name == name) return &g;
            return nullptr;
        }

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
                if (!declared_) fail(tok_, "no qubit register declared");
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
                    fail(tok_, "expected the version statement 'OPENQASM 3.0;' but found " + describe(tok_));
                advance();
                if (tok_.kind != TokenKind::Integer && tok_.kind != TokenKind::Float)
                    fail(tok_, "expected a version number but found " + describe(tok_));
                if (tok_.text != "3.0" && tok_.text != "3")
                    fail(tok_, "unsupported version " + quoted(tok_.text) + "; only 3.0 is supported");
                advance();
                expect(TokenKind::Semicolon, "';'");
            }

            void statement() {
                if (tok_.kind != TokenKind::Identifier)
                    fail(tok_, "expected a statement but found " + describe(tok_));
                if (tok_.text == "OPENQASM")
                    fail(tok_, "the version statement must come first and appear only once");
                if (tok_.text == "include") return include();
                if (tok_.text == "qubit") return declaration();
                if (tok_.text == "gphase") return gphase();
                gate();
            }

            void include() {
                advance();
                const Token file = expect(TokenKind::String, "a file name in double quotes");
                if (file.text != "stdgates.inc")
                    fail(file, "cannot include \"" + std::string(file.text) +
                               "\"; only \"stdgates.inc\" is supported");
                expect(TokenKind::Semicolon, "';'");
                stdgates_ = true;
            }

            void declaration() {
                const Token keyword = advance();
                if (declared_) fail(keyword, "only one qubit register is supported");
                expect(TokenKind::LBracket, "'['");
                const Token size = expect(TokenKind::Integer, "a register size");
                const unsigned long long n = integer(size);
                if (n < 1 || n > static_cast<unsigned long long>(kMaxQubits))
                    fail(size, "register size " + std::string(size.text) + " is outside [1, " +
                               std::to_string(kMaxQubits) + "]");
                expect(TokenKind::RBracket, "']'");
                register_ = expect(TokenKind::Identifier, "a register name").text;
                expect(TokenKind::Semicolon, "';'");
                circuit_.n = static_cast<int>(n);
                declared_ = true;
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
                const GateSpec* spec = find_gate(name.text);
                if (!spec) fail(name, "unknown gate or statement " + quoted(name.text));
                if (!spec->builtin && !stdgates_)
                    fail(name, "gate " + quoted(name.text) + " is not defined; include \"stdgates.inc\" first");
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
                const Token reg = expect(TokenKind::Identifier, "a qubit");
                if (!declared_ || reg.text != register_)
                    fail(reg, "undefined register " + quoted(reg.text));
                expect(TokenKind::LBracket, "'['");
                const Token index = expect(TokenKind::Integer, "a qubit index");
                const unsigned long long i = integer(index);
                if (i >= static_cast<unsigned long long>(circuit_.n))
                    fail(index, "qubit index " + std::string(index.text) + " is out of range for register " +
                                quoted(register_) + " of size " + std::to_string(circuit_.n));
                expect(TokenKind::RBracket, "']'");
                return static_cast<uint8_t>(i);
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
            std::string_view register_;     // name of the qubit register, once declared
            bool declared_ = false;
            bool stdgates_ = false;         // include "stdgates.inc" seen
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
