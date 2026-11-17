#include "lexer.hpp"

#include "qcanopy/parser.hpp"

#include <cstdio>
#include <string>

namespace qcanopy::qasm {

    namespace {
        // Written out rather than taken from <cctype>, whose functions depend on the locale
        // and are undefined for negative char values (any byte >= 0x80).
        bool is_digit(char c) { return c >= '0' && c <= '9'; }
        bool is_ident_start(char c) {
            return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_';
        }
        bool is_ident_char(char c) { return is_ident_start(c) || is_digit(c); }
        bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }

        // Printable ASCII is quoted; anything else is shown as its byte value.
        std::string describe(char c) {
            const unsigned char u = static_cast<unsigned char>(c);
            if (u >= 0x20 && u < 0x7f) return std::string("character '") + c + "'";
            char buf[16];
            std::snprintf(buf, sizeof buf, "byte 0x%02x", static_cast<unsigned>(u));
            return buf;
        }
    }

    Lexer::Lexer(std::string_view source) : src_(source) {}

    char Lexer::peek(std::size_t ahead) const {
        return pos_ + ahead < src_.size() ? src_[pos_ + ahead] : '\0';
    }

    void Lexer::advance() {
        if (src_[pos_] == '\n') {
            ++line_;
            column_ = 1;
        } else {
            ++column_;
        }
        ++pos_;
    }

    Token Lexer::next() {
        while (pos_ < src_.size() && is_space(src_[pos_])) advance();

        Token t{TokenKind::End, src_.substr(pos_, 0), line_, column_};
        if (pos_ == src_.size()) return t;

        const std::size_t start = pos_;
        const char c = src_[pos_];

        if (is_ident_start(c)) {
            while (is_ident_char(peek())) advance();
            t.kind = TokenKind::Identifier;
        } else if (is_digit(c) || (c == '.' && is_digit(peek(1)))) {
            t.kind = number(t);
        } else if (c == '"') {
            advance();
            while (peek() != '"') {
                if (pos_ == src_.size() || peek() == '\n' || peek() == '\r')
                    throw ParseError(t.line, t.column, "unterminated string");
                advance();
            }
            t.kind = TokenKind::String;
            t.text = src_.substr(start + 1, pos_ - start - 1);
            advance();                                  // closing quote
            return t;
        } else {
            switch (c) {
                case '(': t.kind = TokenKind::LParen; break;
                case ')': t.kind = TokenKind::RParen; break;
                case '[': t.kind = TokenKind::LBracket; break;
                case ']': t.kind = TokenKind::RBracket; break;
                case ',': t.kind = TokenKind::Comma; break;
                case ';': t.kind = TokenKind::Semicolon; break;
                case '-': t.kind = TokenKind::Minus; break;
                default:
                    throw ParseError(t.line, t.column, "unexpected " + describe(c));
            }
            advance();
        }

        t.text = src_.substr(start, pos_ - start);
        return t;
    }

    TokenKind Lexer::number(const Token& t) {
        TokenKind kind = TokenKind::Integer;

        while (is_digit(peek())) advance();
        if (peek() == '.') {
            kind = TokenKind::Float;
            advance();
            while (is_digit(peek())) advance();
        }
        if (peek() == 'e' || peek() == 'E') {
            kind = TokenKind::Float;
            advance();
            if (peek() == '+' || peek() == '-') advance();
            if (!is_digit(peek()))
                throw ParseError(t.line, t.column, "exponent has no digits");
            while (is_digit(peek())) advance();
        }

        // "3q" and "1.2.3" are one bad number, not a number followed by something else.
        if (is_ident_char(peek()) || peek() == '.')
            throw ParseError(t.line, t.column, "malformed number");
        return kind;
    }

}
