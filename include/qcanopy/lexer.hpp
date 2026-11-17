// qcanopy/lexer.hpp — first stage of the parser: splits OpenQASM source text into tokens
// (names, numbers, strings, punctuation) so the parser never looks at raw characters.
// To read a circuit, use parse_qasm (qcanopy/parser.hpp) rather than the lexer directly.
//
// 1. Token::text is a view into the source, so the source must outlive its tokens.
// 2. Positions are 1-based. A column is a byte offset within the line, so a tab counts as one.
// 3. Keywords are not separated from names: `qubit`, `rz` and `q` are all Identifier tokens.
// 4. A sign is never part of a number: "-1.5" is Minus followed by Float.
// 5. A number with a point or an exponent is a Float ("5." and ".5" included), anything else
//    an Integer. The value is not converted here.
// 6. Strings are double-quoted and end on the same line. Token::text excludes the quotes.
// 7. Comments are skipped like whitespace: "//" to the end of the line, and "/*" to the
//    next "*/". Block comments do not nest.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace qcanopy::qasm {

    enum class TokenKind : uint8_t {
        Identifier,
        Integer,
        Float,
        String,
        LParen, RParen,         // ( )
        LBracket, RBracket,     // [ ]
        Comma,
        Semicolon,
        Plus, Minus,
        Star, Slash,
        DoubleStar,             // **
        End,                    // end of the source
    };

    struct Token {
        TokenKind kind = TokenKind::End;
        std::string_view text;
        uint32_t line = 1;
        uint32_t column = 1;
    };

    class Lexer {
    public:
        explicit Lexer(std::string_view source);

        // Returns the next token. Once the source is exhausted every call returns End,
        // positioned just past the last character. Throws ParseError on a character or
        // literal that cannot start or continue a token, or a block comment left open.
        Token next();

    private:
        // The character `ahead` positions past the current one, or '\0' beyond the end.
        char peek(std::size_t ahead = 0) const;
        void advance();

        // Skips whitespace and comments up to the next token or the end of the source.
        void skip_blanks();

        // Consumes a number that starts at the current position; t supplies the error position.
        TokenKind number(const Token& t);

        std::string_view src_;
        std::size_t pos_ = 0;
        uint32_t line_ = 1;
        uint32_t column_ = 1;
    };

}
