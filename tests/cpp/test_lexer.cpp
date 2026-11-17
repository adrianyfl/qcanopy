// Tests for the lexer (src/parser/lexer.hpp).
#include "../../src/parser/lexer.hpp"
#include "qcanopy/parser.hpp"

#include <cstdio>
#include <string>
#include <string_view>
#include <vector>

using namespace qcanopy;
using qasm::Lexer;
using qasm::Token;

namespace {
    int failures = 0;

#define CHECK(cond) \
    do { if (!(cond)) { std::printf("%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)

    using K = qasm::TokenKind;

    // Every token up to, but not including, End.
    std::vector<Token> lex(std::string_view src) {
        Lexer lexer(src);
        std::vector<Token> out;
        for (Token t = lexer.next(); t.kind != K::End; t = lexer.next()) out.push_back(t);
        return out;
    }

    struct Expected {
        K kind;
        std::string_view text;
    };

    bool lexes_to(std::string_view src, const std::vector<Expected>& want) {
        const std::vector<Token> got = lex(src);
        if (got.size() != want.size()) return false;
        for (std::size_t i = 0; i < got.size(); ++i)
            if (got[i].kind != want[i].kind || got[i].text != want[i].text) return false;
        return true;
    }

    // True if lexing src throws a ParseError positioned at (line, column).
    bool fails_at(std::string_view src, uint32_t line, uint32_t column) {
        try {
            (void)lex(src);
        } catch (const ParseError& e) {
            return e.line() == line && e.column() == column;
        }
        return false;
    }

    // The statement shapes qcanopy_generate --format qasm writes.
    void test_generator_program() {
        CHECK(lexes_to("OPENQASM 3.0;", {{K::Identifier, "OPENQASM"}, {K::Float, "3.0"}, {K::Semicolon, ";"}}));
        CHECK(lexes_to("include \"stdgates.inc\";",
                       {{K::Identifier, "include"}, {K::String, "stdgates.inc"}, {K::Semicolon, ";"}}));
        CHECK(lexes_to("qubit[4] q;",
                       {{K::Identifier, "qubit"}, {K::LBracket, "["}, {K::Integer, "4"}, {K::RBracket, "]"},
                        {K::Identifier, "q"}, {K::Semicolon, ";"}}));
        CHECK(lexes_to("gphase(4.3196898986859651);",
                       {{K::Identifier, "gphase"}, {K::LParen, "("}, {K::Float, "4.3196898986859651"},
                        {K::RParen, ")"}, {K::Semicolon, ";"}}));
        CHECK(lexes_to("rz(-0.78539816339744828) q[0];",
                       {{K::Identifier, "rz"}, {K::LParen, "("}, {K::Minus, "-"},
                        {K::Float, "0.78539816339744828"}, {K::RParen, ")"}, {K::Identifier, "q"},
                        {K::LBracket, "["}, {K::Integer, "0"}, {K::RBracket, "]"}, {K::Semicolon, ";"}}));
        CHECK(lexes_to("U(0.25,0.5,0.75) q[3];",
                       {{K::Identifier, "U"}, {K::LParen, "("}, {K::Float, "0.25"}, {K::Comma, ","},
                        {K::Float, "0.5"}, {K::Comma, ","}, {K::Float, "0.75"}, {K::RParen, ")"},
                        {K::Identifier, "q"}, {K::LBracket, "["}, {K::Integer, "3"}, {K::RBracket, "]"},
                        {K::Semicolon, ";"}}));
        CHECK(lexes_to("cx q[1],q[0];",
                       {{K::Identifier, "cx"}, {K::Identifier, "q"}, {K::LBracket, "["}, {K::Integer, "1"},
                        {K::RBracket, "]"}, {K::Comma, ","}, {K::Identifier, "q"}, {K::LBracket, "["},
                        {K::Integer, "0"}, {K::RBracket, "]"}, {K::Semicolon, ";"}}));
    }

    void test_identifiers() {
        CHECK(lexes_to("sdg _tmp q0 CX", {{K::Identifier, "sdg"}, {K::Identifier, "_tmp"},
                                         {K::Identifier, "q0"}, {K::Identifier, "CX"}}));
    }

    void test_numbers() {
        CHECK(lexes_to("0", {{K::Integer, "0"}}));
        CHECK(lexes_to("42", {{K::Integer, "42"}}));
        CHECK(lexes_to("3.0", {{K::Float, "3.0"}}));
        CHECK(lexes_to("5.", {{K::Float, "5."}}));
        CHECK(lexes_to(".5", {{K::Float, ".5"}}));
        CHECK(lexes_to("1e-05", {{K::Float, "1e-05"}}));
        CHECK(lexes_to("1.5E+20", {{K::Float, "1.5E+20"}}));
        CHECK(lexes_to("2e3", {{K::Float, "2e3"}}));
        // The sign is its own token.
        CHECK(lexes_to("-1.5", {{K::Minus, "-"}, {K::Float, "1.5"}}));
        CHECK(lexes_to("--2", {{K::Minus, "-"}, {K::Minus, "-"}, {K::Integer, "2"}}));
    }

    void test_strings() {
        CHECK(lexes_to("\"a b.inc\"", {{K::String, "a b.inc"}}));
        CHECK(lexes_to("\"\"", {{K::String, ""}}));
        // Nothing inside a string is tokenized.
        CHECK(lexes_to("\"q[0]; // x\" ;", {{K::String, "q[0]; // x"}, {K::Semicolon, ";"}}));
    }

    void test_positions() {
        // CRLF line ending, then a tab before the second statement.
        const std::vector<Token> t = lex("h q[0];\r\n\tcx q[1],q[0];");
        CHECK(t.size() == 17);
        if (t.size() != 17) return;
        CHECK(t[0].line == 1 && t[0].column == 1);      // h
        CHECK(t[1].line == 1 && t[1].column == 3);      // q
        CHECK(t[5].line == 1 && t[5].column == 7);      // ;
        CHECK(t[6].line == 2 && t[6].column == 2);      // cx
        CHECK(t[11].line == 2 && t[11].column == 9);    // ,
        CHECK(t[16].line == 2 && t[16].column == 14);   // ;

        // A string token is positioned at its opening quote.
        const std::vector<Token> s = lex("include \"x.inc\";");
        CHECK(s.size() == 3 && s[1].column == 9 && s[2].column == 16);
    }

    void test_end() {
        Lexer empty("");
        CHECK(empty.next().kind == K::End);
        CHECK(empty.next().kind == K::End);

        Lexer blank(" \t\r\n\n  ");
        const Token e = blank.next();
        CHECK(e.kind == K::End && e.text.empty());
        CHECK(e.line == 3 && e.column == 3);

        Lexer one("h\n");
        CHECK(one.next().kind == K::Identifier);
        const Token end = one.next();
        CHECK(end.kind == K::End && end.line == 2 && end.column == 1);
        CHECK(one.next().kind == K::End);
    }

    void test_errors() {
        CHECK(fails_at("h q[0];\n  @", 2, 3));
        CHECK(fails_at("// comment", 1, 1));
        CHECK(fails_at("rz(pi/2) q[0];", 1, 6));           // '/': no expressions
        CHECK(fails_at("rz(+1) q[0];", 1, 4));
        CHECK(fails_at("x q[0]; \xc3\xa9", 1, 9));         // UTF-8 bytes
        CHECK(fails_at(std::string_view("h\0q", 3), 1, 2));

        // Bad numbers are reported where the number starts.
        CHECK(fails_at("rz(1e) q[0];", 1, 4));
        CHECK(fails_at("rz(1e+) q[0];", 1, 4));
        CHECK(fails_at("rz(1.2.3) q[0];", 1, 4));
        CHECK(fails_at("q[3q]", 1, 3));
        CHECK(fails_at("0x1F", 1, 1));
        CHECK(fails_at(".", 1, 1));

        // A string must close on its own line.
        CHECK(fails_at("include \"stdgates.inc;", 1, 9));
        CHECK(fails_at("include \"stdgates\n.inc\";", 1, 9));

        try {
            (void)lex("h q[0];\n  @");
            CHECK(false);
        } catch (const ParseError& e) {
            CHECK(std::string(e.what()) == "2:3: unexpected character '@'");
        }
        try {
            (void)lex("\x01");
            CHECK(false);
        } catch (const std::invalid_argument& e) {      // ParseError is an invalid_argument
            CHECK(std::string(e.what()) == "1:1: unexpected byte 0x01");
        }
    }
} // namespace

int main() {
    test_generator_program();
    test_identifiers();
    test_numbers();
    test_strings();
    test_positions();
    test_end();
    test_errors();
    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("all lexer tests passed\n");
    return 0;
}
