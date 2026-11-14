// Tests for normalize(): canonical forms, fusion, cancellation, and exact unitary preservation.
#include "qcanopy/dag.hpp"
#include "qcanopy/normalize.hpp"

#include <cmath>
#include <cstdio>
#include <random>

using namespace qcanopy;

namespace {
    int failures = 0;

#define CHECK(cond) \
    do { if (!(cond)) { std::printf("%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)

    constexpr double kPi = 3.14159265358979323846;

    // Dense unitary of c (column-major: u[col * dim + row]), including global phase.
    std::vector<cplx> unitary(const Circuit& c) {
        const std::size_t dim = std::size_t{1} << c.n;
        std::vector<cplx> u(dim * dim);
        for (std::size_t col = 0; col < dim; ++col) {
            cplx* v = &u[col * dim];
            v[col] = std::polar(1.0, c.global_phase);
            for (const Gate& g : c.gates) {
                if (g.nq == 1) {
                    const std::size_t b = std::size_t{1} << g.q[0];
                    for (std::size_t i = 0; i < dim; ++i) {
                        if (i & b) continue;
                        const cplx a0 = v[i], a1 = v[i | b];
                        v[i]     = g.at(0, 0) * a0 + g.at(0, 1) * a1;
                        v[i | b] = g.at(1, 0) * a0 + g.at(1, 1) * a1;
                    }
                } else {
                    const std::size_t b0 = std::size_t{1} << g.q[0], b1 = std::size_t{1} << g.q[1];
                    for (std::size_t i = 0; i < dim; ++i) {
                        if (i & (b0 | b1)) continue;
                        const std::size_t idx[4] = {i, i | b0, i | b1, i | b0 | b1};   // local index (b_q1 << 1) | b_q0
                        cplx old[4], out[4] = {};
                        for (int k = 0; k < 4; ++k) old[k] = v[idx[k]];
                        for (int r = 0; r < 4; ++r)
                            for (int k = 0; k < 4; ++k) out[r] += g.at(r, k) * old[k];
                        for (int k = 0; k < 4; ++k) v[idx[k]] = out[k];
                    }
                }
            }
        }
        return u;
    }

    bool same_unitary(const Circuit& a, const Circuit& b, double tol = 1e-9) {
        const auto ua = unitary(a), ub = unitary(b);
        for (std::size_t i = 0; i < ua.size(); ++i)
            if (std::abs(ua[i] - ub[i]) > tol) return false;
        return true;
    }

    bool same_circuit(const Circuit& a, const Circuit& b) {
        if (a.n != b.n || a.size() != b.size()) return false;
        if (std::abs(std::remainder(a.global_phase - b.global_phase, 2 * kPi)) > 1e-9) return false;
        for (std::size_t i = 0; i < a.size(); ++i)
            if (!same_gate(a.gates[i], b.gates[i], 1e-9)) return false;
        return true;
    }

    // Normalizes c and checks the unitary is unchanged.
    Circuit normalized(const Circuit& c) {
        CircuitDAG dag(c);
        normalize(dag);
        Circuit out = dag.to_circuit();
        CHECK(same_unitary(c, out));
        return out;
    }

    Gate rz(uint8_t q, double t) { return make_gate(GateKind::RZ, q, {t, 0, 0}); }
    Gate rx(uint8_t q, double t) { return make_gate(GateKind::RX, q, {t, 0, 0}); }

    void test_to_circuit_roundtrip() {
        Circuit c(3);
        c.add(make_gate(GateKind::H, 0)).add(make_gate(GateKind::CX, 0, 1)).add(make_gate(GateKind::T, 2))
         .add(make_gate(GateKind::CX, 1, 2)).add(make_gate(GateKind::H, 0));
        CHECK(same_unitary(c, CircuitDAG(c).to_circuit()));
        CHECK(CircuitDAG(c).to_circuit().size() == c.size());
    }

    void test_z_family() {
        Circuit c(1);
        c.add(make_gate(GateKind::T, 0)).add(make_gate(GateKind::T, 0));
        Circuit out = normalized(c);
        CHECK(out.size() == 1 && out.gates[0].kind == GateKind::S);

        Circuit c2(1);
        c2.add(rz(0, kPi / 4));                    // RZ(pi/4) = e^{-i pi/8} T
        out = normalized(c2);
        CHECK(out.size() == 1 && out.gates[0].kind == GateKind::T);
        CHECK(std::abs(out.global_phase + kPi / 8) < 1e-12);

        Circuit c3(1);
        c3.add(rz(0, 0.3)).add(rz(0, 0.4));
        out = normalized(c3);
        CHECK(out.size() == 1 && out.gates[0].kind == GateKind::P && std::abs(out.gates[0].params[0] - 0.7) < 1e-12);

        Circuit c4(1);
        c4.add(make_gate(GateKind::T, 0)).add(make_gate(GateKind::Tdg, 0)).add(make_gate(GateKind::I, 0));
        CHECK(normalized(c4).size() == 0);

        Circuit c5(1);
        c5.add(make_gate(GateKind::U, 0, {0, 0.2, 0.5}));
        out = normalized(c5);
        CHECK(out.size() == 1 && out.gates[0].kind == GateKind::P);
    }

    void test_x_family() {
        Circuit c(1);
        c.add(make_gate(GateKind::SX, 0)).add(make_gate(GateKind::SX, 0));
        Circuit out = normalized(c);
        CHECK(out.size() == 1 && out.gates[0].kind == GateKind::X);

        Circuit c2(1);
        c2.add(rx(0, 3 * kPi));                    // RX(3pi) = -RX(pi) = i X
        out = normalized(c2);
        CHECK(out.size() == 1 && out.gates[0].kind == GateKind::X);

        Circuit c3(1);
        c3.add(rx(0, 2 * kPi));                    // -I
        out = normalized(c3);
        CHECK(out.size() == 0 && std::abs(std::abs(out.global_phase) - kPi) < 1e-12);
    }

    void test_cancellation() {
        Circuit c(2);
        c.add(make_gate(GateKind::H, 0)).add(make_gate(GateKind::CX, 0, 1)).add(rz(1, 0.3))
         .add(rz(1, -0.3)).add(make_gate(GateKind::CX, 0, 1)).add(make_gate(GateKind::H, 0));
        CHECK(normalized(c).size() == 0);          // needs fuse -> cancel -> cancel to cascade

        Circuit c2(2);
        c2.add(make_gate(GateKind::CZ, 1, 0)).add(make_gate(GateKind::CZ, 0, 1));
        CHECK(normalized(c2).size() == 0);

        Circuit c3(2);                             // CX(0,1) CX(1,0) is not an inverse pair
        c3.add(make_gate(GateKind::CX, 0, 1)).add(make_gate(GateKind::CX, 1, 0));
        CHECK(normalized(c3).size() == 2);

        Circuit c4(2);                             // the CX blocks fusion on wire 0
        c4.add(rz(0, 0.3)).add(make_gate(GateKind::CX, 0, 1)).add(rz(0, 0.4));
        CHECK(normalized(c4).size() == 3);
    }

    void test_symmetric_order() {
        Circuit c(3);
        c.add(make_gate(GateKind::CZ, 2, 0)).add(make_gate(GateKind::SWAP, 1, 0));
        Circuit out = normalized(c);
        for (const Gate& g : out.gates) CHECK(g.q[0] < g.q[1]);
    }

    // Two spellings of the same circuit normalize to the same gate list.
    void test_canonical_spellings() {
        Circuit a(3), b(3);
        a.add(make_gate(GateKind::H, 0)).add(make_gate(GateKind::T, 1)).add(make_gate(GateKind::CZ, 0, 1))
         .add(make_gate(GateKind::SX, 2)).add(make_gate(GateKind::SX, 2)).add(make_gate(GateKind::CX, 1, 2));
        b.add(rz(1, kPi / 8)).add(make_gate(GateKind::P, 1, {kPi / 8, 0, 0})).add(rx(2, kPi))
         .add(make_gate(GateKind::H, 0)).add(make_gate(GateKind::CZ, 1, 0)).add(make_gate(GateKind::CX, 1, 2));
        b.global_phase = kPi / 16 + kPi / 2;       // undo RZ's e^{-i pi/16} and RX(pi) = -i X
        CHECK(same_unitary(a, b));
        CHECK(same_circuit(normalized(a), normalized(b)));
    }

    // Random circuits with planted inverses: unitary preserved, result idempotent.
    void test_random() {
        std::mt19937 rng(12345);
        const GateKind kinds1[] = {
            GateKind::I, GateKind::H, GateKind::X, GateKind::Y, GateKind::Z, GateKind::S, GateKind::Sdg,
            GateKind::T, GateKind::Tdg, GateKind::SX, GateKind::SXdg, GateKind::RX, GateKind::RY,
            GateKind::RZ, GateKind::P, GateKind::U,
        };
        const GateKind kinds2[] = {GateKind::CX, GateKind::CZ, GateKind::SWAP};
        auto angle = [&] {
            const int k = static_cast<int>(rng() % 17) - 8;
            return (rng() % 2) ? k * kPi / 4 : std::uniform_real_distribution<double>(-7, 7)(rng);
        };

        for (int trial = 0; trial < 300; ++trial) {
            Circuit c(3);
            for (int i = 0; i < 40; ++i) {
                const uint8_t q0 = static_cast<uint8_t>(rng() % 3);
                const uint8_t q1 = static_cast<uint8_t>((q0 + 1 + rng() % 2) % 3);
                if (!c.gates.empty() && rng() % 5 == 0) { c.add(inverse(c.gates.back())); continue; }
                if (rng() % 4 == 0) c.add(make_gate(kinds2[rng() % 3], q0, q1));
                else c.add(make_gate(kinds1[rng() % 16], q0, {angle(), angle(), angle()}));
            }
            const Circuit once = normalized(c);
            CHECK(once.size() <= c.size());
            CHECK(same_circuit(once, normalized(once)));
        }
    }
}

int main() {
    test_to_circuit_roundtrip();
    test_z_family();
    test_x_family();
    test_cancellation();
    test_symmetric_order();
    test_canonical_spellings();
    test_random();
    if (failures) std::printf("%d check(s) failed\n", failures);
    else std::printf("all tests passed\n");
    return failures ? 1 : 0;
}
