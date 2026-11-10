// Private strict JSON support. No network, generated code, or third-party dependency.
#pragma once
#include <cmath>
#include <iomanip>
#include <istream>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qcanopy::testgen::json {
struct Value {
    enum Type { Null, Boolean, Number, String, Array, Object } type = Null;
    double number = 0;
    std::string string;
    std::vector<Value> array;
    std::map<std::string,Value> object;
    const Value& field(const std::string& key) const {
        if (type != Object || !object.count(key)) throw std::invalid_argument("missing JSON field: " + key);
        return object.at(key);
    }
    void keys(std::initializer_list<const char*> expected) const {
        if (type != Object || object.size() != expected.size()) throw std::invalid_argument("unexpected JSON object fields");
        for (auto k : expected) (void)field(k);
    }
    double real() const {
        if (type != Number || !std::isfinite(number)) throw std::invalid_argument("expected finite JSON number");
        return number;
    }
    int integer(int max) const {
        const double x = real();
        if (x < 0 || x > max || std::floor(x) != x) throw std::invalid_argument("JSON integer out of range");
        return static_cast<int>(x);
    }
    const std::vector<Value>& items() const {
        if (type != Array) throw std::invalid_argument("expected JSON array");
        return array;
    }
    const std::string& str() const {
        if (type != String) throw std::invalid_argument("expected JSON string");
        return string;
    }
};
class Parser {
    std::istream& in;
    [[noreturn]] void fail() { throw std::invalid_argument("malformed JSON circuit"); }
    void ws() { while (in.peek() == ' ' || in.peek() == '\n' || in.peek() == '\r' || in.peek() == '\t') in.get(); }
    void take(char c) { if (in.get() != c) fail(); }
    static bool digit(int c) { return c >= '0' && c <= '9'; }
    unsigned hex4() {
        unsigned v = 0;
        for (int i = 0; i < 4; ++i) {
            const int c = in.get(); unsigned n;
            if (c >= '0' && c <= '9') n = static_cast<unsigned>(c-'0');
            else if (c >= 'a' && c <= 'f') n = static_cast<unsigned>(c-'a'+10);
            else if (c >= 'A' && c <= 'F') n = static_cast<unsigned>(c-'A'+10);
            else { fail(); }
            v = v * 16 + n;
        }
        return v;
    }
    static void utf8(std::string& s, unsigned cp) {
        if (cp < 0x80) s += static_cast<char>(cp);
        else if (cp < 0x800) { s += static_cast<char>(0xc0|(cp>>6)); s += static_cast<char>(0x80|(cp&63)); }
        else if (cp < 0x10000) {
            s += static_cast<char>(0xe0|(cp>>12)); s += static_cast<char>(0x80|((cp>>6)&63)); s += static_cast<char>(0x80|(cp&63));
        } else {
            s += static_cast<char>(0xf0|(cp>>18)); s += static_cast<char>(0x80|((cp>>12)&63));
            s += static_cast<char>(0x80|((cp>>6)&63)); s += static_cast<char>(0x80|(cp&63));
        }
    }
    std::string string() {
        take('"'); std::string out;
        for (;;) {
            const int c = in.get();
            if (c == '"') return out;
            if (c < 32) fail();
            if (c != '\\') { out += static_cast<char>(c); continue; }
            switch (in.get()) {
                case '"': out += '"'; break; case '\\': out += '\\'; break; case '/': out += '/'; break;
                case 'b': out += '\b'; break; case 'f': out += '\f'; break; case 'n': out += '\n'; break;
                case 'r': out += '\r'; break; case 't': out += '\t'; break;
                case 'u': {
                    unsigned cp = hex4();
                    if (cp >= 0xd800 && cp <= 0xdbff) {
                        take('\\'); take('u'); const unsigned low = hex4();
                        if (low < 0xdc00 || low > 0xdfff) fail();
                        cp = 0x10000 + ((cp-0xd800)<<10) + low-0xdc00;
                    } else if (cp >= 0xdc00 && cp <= 0xdfff) fail();
                    utf8(out,cp); break;
                }
                default: fail();
            }
        }
    }
    Value value(int depth) {
        if (depth > 64) fail();
        ws(); Value v;
        if (in.peek() == '{') {
            in.get(); v.type = Value::Object; ws();
            if (in.peek() == '}') { in.get(); return v; }
            for (;;) {
                ws(); auto key = string(); ws(); take(':');
                auto child = value(depth+1);
                if (!v.object.emplace(key,std::move(child)).second) fail();
                ws(); const int c = in.get(); if (c == '}') return v; if (c != ',') fail();
            }
        }
        if (in.peek() == '[') {
            in.get(); v.type = Value::Array; ws();
            if (in.peek() == ']') { in.get(); return v; }
            for (;;) {
                v.array.push_back(value(depth+1)); ws();
                const int c = in.get(); if (c == ']') return v; if (c != ',') fail();
            }
        }
        if (in.peek() == '"') { v.type = Value::String; v.string = string(); return v; }
        if (in.peek() == 't' || in.peek() == 'f' || in.peek() == 'n') {
            const char* word = in.peek() == 't' ? "true" : in.peek() == 'f' ? "false" : "null";
            v.type = *word == 'n' ? Value::Null : Value::Boolean;
            while (*word) take(*word++);
            return v;
        }
        std::string token;
        auto consume = [&] { token += static_cast<char>(in.get()); };
        if (in.peek() == '-') consume();
        if (in.peek() == '0') consume();
        else {
            if (!digit(in.peek())) fail();
            while (digit(in.peek())) consume();
        }
        if (in.peek() == '.') {
            consume(); if (!digit(in.peek())) fail(); while (digit(in.peek())) consume();
        }
        if (in.peek() == 'e' || in.peek() == 'E') {
            consume(); if (in.peek() == '+' || in.peek() == '-') consume();
            if (!digit(in.peek())) fail();
            while (digit(in.peek())) consume();
        }
        std::istringstream num(token); num.imbue(std::locale::classic());
        v.type = Value::Number; num >> v.number;
        if (!num || !std::isfinite(v.number)) fail();
        return v;
    }
public:
    explicit Parser(std::istream& stream) : in(stream) {}
    Value parse() {
        Value v = value(0); ws();
        if (in.peek() != std::char_traits<char>::eof() || in.bad()) fail();
        return v;
    }
};
inline std::string quote(const std::string& s) {
    std::ostringstream out; out << '"';
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<unsigned>(c) << std::dec;
        else out << c;
    }
    out << '"'; return out.str();
}
inline void numbers(std::ostream& s) { s.imbue(std::locale::classic()); s << std::setprecision(std::numeric_limits<double>::max_digits10); }
} // namespace qcanopy::testgen::json
