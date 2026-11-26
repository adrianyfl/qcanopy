// Tests for EquivalenceChecker, miter(), the Choi-state checker, and Simulator::execute_each
// (streaming + early stop).
#include "qcanopy/equiv.hpp"
#include "qcanopy/stimulus.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

using namespace qcanopy;

namespace {
    int failures = 0;

#define CHECK(cond) \
    do { if (!(cond)) { std::printf("%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)

    using K = GateKind;

    bool same_state(const StateVector& a, const StateVector& b, double tol) {
        if (a.amplitudes.size() != b.amplitudes.size()) return false;
        for (std::size_t i = 0; i < a.amplitudes.size(); ++i)
            if (std::abs(a.amplitudes[i] - b.amplitudes[i]) > tol) return false;
        return true;
    }

    template <class F>
    bool throws_invalid(F f) {
        try { f(); } catch (const std::invalid_argument&) { return true; }
        return false;
    }

    Stimulus mixed(int n, uint64_t seed) {
        Stimulus s = haar_random_states(n, 4, seed);
        s.append(random_product_states(n, 4, seed)).append(all_basis_states(n));
        return s;
    }

    // A fixed 3-qubit circuit with a mix of fixed, parameterized, and 2-qubit gates.
    Circuit sample_circuit() {
        Circuit c(3);
        c.add(make_gate(K::H, 0)).add(make_gate(K::CX, 0, 1)).add(make_gate(K::RY, 2, {0.37, 0, 0}))
         .add(make_gate(K::T, 1)).add(make_gate(K::CZ, 1, 2)).add(make_gate(K::U, 0, {0.3, -1.1, 2.4}))
         .add(make_gate(K::SWAP, 0, 2)).add(make_gate(K::SXdg, 1)).add(make_gate(K::RZ, 2, {-0.8, 0, 0}));
        return c;
    }

    void test_execute_each() {
        Simulator sim;
        const Circuit c = sample_circuit();
        const Stimulus stim = haar_random_states(3, 5, 1);

        // No stop: every index visited once, outputs match single-state execution.
        std::vector<std::size_t> seen;
        sim.execute_each(c, stim, [&](std::size_t i, const StateVector& out) {
            seen.push_back(i);
            CHECK(same_state(out, sim.execute(c, stim.states[i]), 0.0));
            return Control::Continue;
        });
        CHECK((seen == std::vector<std::size_t>{0, 1, 2, 3, 4}));

        // Stop at index 2 cancels everything after it.
        seen.clear();
        sim.execute_each(c, stim, [&](std::size_t i, const StateVector&) {
            seen.push_back(i);
            return i == 2 ? Control::Stop : Control::Continue;
        });
        CHECK((seen == std::vector<std::size_t>{0, 1, 2}));

        CHECK(throws_invalid([&] { sim.execute_each(c, haar_random_states(2, 1, 0), [](std::size_t, const StateVector&) { return Control::Continue; }); }));
    }

    void test_miter() {
        Simulator sim;
        const Circuit a = sample_circuit();
        // A^-1 A = I: every input comes back unchanged.
        const Circuit m = miter(a, a);
        CHECK(m.size() == 2 * a.size());
        for (const auto& psi : haar_random_states(3, 3, 9).states) CHECK(same_state(sim.execute(m, psi), psi, 1e-12));
        CHECK(throws_invalid([&] { (void)miter(a, Circuit(2)); }));
    }

    void test_equivalent_identities() {
        const EquivalenceChecker checker;

        // H X H = Z
        {
            Circuit a(1), b(1);
            a.add(make_gate(K::H, 0)).add(make_gate(K::X, 0)).add(make_gate(K::H, 0));
            b.add(make_gate(K::Z, 0));
            const EquivResult r = checker.check(a, b, mixed(1, 1));
            CHECK(r.equivalent);
            CHECK(r.mismatch == Mismatch::None);
            CHECK(r.stimuli_checked == 4 + 4 + 2);
            CHECK(r.deviation <= checker.tolerance());
            CHECK(std::abs(r.global_phase - cplx{1.0, 0.0}) < 1e-12);
        }
        // CX(0,1) = (H x H) CX(1,0) (H x H)
        {
            Circuit a(2), b(2);
            a.add(make_gate(K::CX, 0, 1));
            b.add(make_gate(K::H, 0)).add(make_gate(K::H, 1)).add(make_gate(K::CX, 1, 0))
             .add(make_gate(K::H, 0)).add(make_gate(K::H, 1));
            CHECK(checker.check(a, b, mixed(2, 2)).equivalent);
        }
        // SWAP = CX CX CX
        {
            Circuit a(2), b(2);
            a.add(make_gate(K::SWAP, 0, 1));
            b.add(make_gate(K::CX, 0, 1)).add(make_gate(K::CX, 1, 0)).add(make_gate(K::CX, 0, 1));
            CHECK(checker.check(a, b, mixed(2, 3)).equivalent);
        }
        // Redundant pairs inserted into a larger circuit.
        {
            const Circuit a = sample_circuit();
            Circuit b(3);
            for (std::size_t g = 0; g < a.size(); ++g) {
                b.add(a.gates[g]);
                if (g == 2) b.add(make_gate(K::H, 1)).add(make_gate(K::H, 1));
                if (g == 5) b.add(make_gate(K::CX, 2, 0)).add(make_gate(K::CX, 2, 0));
            }
            CHECK(checker.check(a, b, mixed(3, 4)).equivalent);
        }
    }

    void test_global_phase() {
        const EquivalenceChecker checker;

        // RZ(t) = e^{-it/2} P(t): equivalent, and the phase is reported.
        const double t = 0.9;
        Circuit a(1), b(1);
        a.add(make_gate(K::RZ, 0, {t, 0, 0}));
        b.add(make_gate(K::P, 0, {t, 0, 0}));
        const EquivResult r = checker.check(a, b, mixed(1, 5));
        CHECK(r.equivalent);
        CHECK(std::abs(r.global_phase - std::polar(1.0, -t / 2)) < 1e-12);

        // The circuit-level global_phase field is honored too.
        Circuit c = sample_circuit();
        c.global_phase = 0.7;
        const EquivResult r2 = checker.check(c, sample_circuit(), mixed(3, 6));
        CHECK(r2.equivalent);
        CHECK(std::abs(r2.global_phase - std::polar(1.0, 0.7)) < 1e-12);
    }

    void test_not_equivalent() {
        const EquivalenceChecker checker;

        // CX vs CZ: the first Haar stimulus catches it, and the run stops there.
        {
            Circuit a(2), b(2);
            a.add(make_gate(K::CX, 0, 1));
            b.add(make_gate(K::CZ, 0, 1));
            const EquivResult r = checker.check(a, b, mixed(2, 7));
            CHECK(!r.equivalent);
            CHECK(r.mismatch == Mismatch::State);
            CHECK(r.failing_index == 0);
            CHECK(r.failing_kind == StimulusKind::Haar);
            CHECK(r.stimuli_checked == 1);
            CHECK(r.deviation > checker.tolerance());
        }
        // A tiny angle error is still caught.
        {
            Circuit a = sample_circuit(), b = sample_circuit();
            b.gates[2] = make_gate(K::RY, 2, {0.37 + 1e-6, 0, 0});
            const EquivResult r = checker.check(a, b, haar_random_states(3, 4, 8));
            CHECK(!r.equivalent);
            CHECK(r.failing_index == 0);
        }
    }

    // Z vs I agree on each basis input up to phase; only the common-phase test separates them.
    void test_phase_mismatch_and_lowest_index() {
        const EquivalenceChecker checker;
        Circuit z(1), id(1);
        z.add(make_gate(K::Z, 0));

        const EquivResult r = checker.check(z, id, all_basis_states(1));
        CHECK(!r.equivalent);
        CHECK(r.mismatch == Mismatch::Phase);
        CHECK(r.failing_index == 1);
        CHECK(r.failing_kind == StimulusKind::Basis);
        CHECK(std::abs(r.deviation - 2.0) < 1e-12);   // |(-1) - (+1)|

        // Stimuli [|0>, |0>, |1>, |1>]: indices 0-1 pass, 2 fails, 3 is never run.
        Stimulus stim(1);
        stim.add(basis_state(1, 0)).add(basis_state(1, 0)).add(basis_state(1, 1)).add(basis_state(1, 1));
        const EquivResult r2 = checker.check(z, id, stim);
        CHECK(!r2.equivalent);
        CHECK(r2.failing_index == 2);
        CHECK(r2.failing_kind == StimulusKind::Custom);
        CHECK(r2.stimuli_checked == 3);
    }

    void test_invalid_inputs() {
        const EquivalenceChecker checker;
        const Circuit a(2), b(3);
        CHECK(throws_invalid([&] { (void)checker.check(a, b, all_basis_states(2)); }));
        CHECK(throws_invalid([&] { (void)checker.check(a, a, all_basis_states(3)); }));
        CHECK(throws_invalid([&] { (void)checker.check(a, a, Stimulus(2)); }));
        CHECK(throws_invalid([] { EquivalenceChecker bad(0.0); }));
        CHECK(throws_invalid([] { EquivalenceChecker bad(-1e-9); }));
    }
    // ---- Choi-state checking ----

    // Every gate kind, including Y and generic U1Q/U2Q, on a 3-qubit register.
    Circuit all_kinds_circuit() {
        const double r = 1.0 / std::sqrt(2.0);
        const cplx i{0.0, 1.0};
        std::array<cplx, 16> u2{};   // a non-real 2-qubit unitary: diag(1, i, -1, e^{i 0.3}) then permuted
        u2[0 * 4 + 1] = 1.0; u2[1 * 4 + 0] = i; u2[2 * 4 + 3] = -1.0; u2[3 * 4 + 2] = std::polar(1.0, 0.3);
        Circuit c(3);
        c.add(make_gate(K::I, 0)).add(make_gate(K::H, 1)).add(make_gate(K::X, 2)).add(make_gate(K::Y, 0))
         .add(make_gate(K::Z, 1)).add(make_gate(K::S, 2)).add(make_gate(K::Sdg, 0)).add(make_gate(K::T, 1))
         .add(make_gate(K::Tdg, 2)).add(make_gate(K::SX, 0)).add(make_gate(K::SXdg, 1))
         .add(make_gate(K::RX, 2, {0.41, 0, 0})).add(make_gate(K::RY, 0, {-1.3, 0, 0}))
         .add(make_gate(K::RZ, 1, {2.2, 0, 0})).add(make_gate(K::P, 2, {0.77, 0, 0}))
         .add(make_gate(K::U, 0, {0.3, -1.1, 2.4})).add(make_gate(K::CX, 0, 2)).add(make_gate(K::CZ, 2, 1))
         .add(make_gate(K::SWAP, 1, 0))
         .add(make_unitary_1q(2, {r, i * r, i * r, r}))
         .add(make_unitary_2q(1, 2, u2));
        return c;
    }

    void test_conjugate() {
        for (const Gate& g : all_kinds_circuit().gates) {
            const Gate c = conjugate(g);
            CHECK(c.nq == g.nq);
            CHECK(c.q == g.q);
            const int d = g.nq == 1 ? 4 : 16;
            for (int k = 0; k < d; ++k) CHECK(std::abs(c.m[k] - std::conj(g.m[k])) < 1e-12);
        }
        Circuit c = sample_circuit();
        c.global_phase = 0.4;
        const Circuit cc = conjugate(c);
        CHECK(cc.size() == c.size());
        CHECK(cc.global_phase == -0.4);
    }

    void test_choi_state() {
        Simulator sim;
        for (int n = 0; n <= 4; ++n) {
            const StateVector phi = choi_state(n);
            CHECK(phi.num_qbits == 2 * n);
            // Matches the H + CX preparation circuit run from |0...0>.
            CHECK(same_state(phi, sim.execute(choi_prep_circuit(n)), 1e-12));
            double norm = 0.0;
            for (const cplx& a : phi.amplitudes) norm += std::norm(a);
            CHECK(std::abs(norm - 1.0) < 1e-12);
        }
        // Qubit k is paired with qubit k + n: amplitude at i | (i << n).
        const StateVector phi = choi_state(2);
        for (std::size_t idx = 0; idx < 16; ++idx) {
            const bool on_diag = (idx & 3) == (idx >> 2);
            CHECK(std::abs(phi.amplitudes[idx] - cplx{on_diag ? 0.5 : 0.0, 0.0}) < 1e-12);
        }
        CHECK(throws_invalid([] { (void)choi_state(-1); }));
        CHECK(throws_invalid([] { (void)choi_state(kMaxQubits / 2 + 1); }));
        CHECK(throws_invalid([] { (void)choi_prep_circuit(kMaxQubits / 2 + 1); }));
    }

    void test_choi_circuit() {
        Simulator sim;
        const Circuit r = all_kinds_circuit();
        const Circuit cc = choi_circuit(r, r);
        CHECK(cc.n == 6);
        CHECK(cc.size() == 2 * r.size());
        // Second-half gates are shifted by n.
        CHECK(cc.gates[r.size()].q[0] == r.gates[0].q[0] + 3);
        // Transpose trick: (R^* (x) R)|Phi> = (I (x) R R^dagger)|Phi> = |Phi>.
        CHECK(same_state(sim.execute(cc, choi_state(3)), choi_state(3), 1e-12));
        // ... and the same with R^* and C applied to a simulator-prepared Choi state.
        Circuit prep = choi_prep_circuit(3);
        for (const Gate& g : cc.gates) prep.add(g);
        CHECK(same_state(sim.execute(prep), choi_state(3), 1e-12));

        CHECK(throws_invalid([] { (void)choi_circuit(Circuit(2), Circuit(3)); }));
        CHECK(throws_invalid([] { (void)choi_circuit(Circuit(17), Circuit(17)); }));
    }

    void test_check_choi_equivalent() {
        const EquivalenceChecker checker;
        // H X H = Z
        {
            Circuit a(1), b(1);
            a.add(make_gate(K::H, 0)).add(make_gate(K::X, 0)).add(make_gate(K::H, 0));
            b.add(make_gate(K::Z, 0));
            const EquivResult r = checker.check_choi(a, b);
            CHECK(r.equivalent);
            CHECK(r.mismatch == Mismatch::None);
            CHECK(r.stimuli_checked == 1);
            CHECK(r.deviation <= checker.tolerance());
            CHECK(std::abs(r.global_phase - cplx{1.0, 0.0}) < 1e-12);
        }
        // SWAP = CX CX CX
        {
            Circuit a(2), b(2);
            a.add(make_gate(K::SWAP, 0, 1));
            b.add(make_gate(K::CX, 0, 1)).add(make_gate(K::CX, 1, 0)).add(make_gate(K::CX, 0, 1));
            CHECK(checker.check_choi(a, b).equivalent);
        }
        // Every gate kind against itself with redundant pairs inserted.
        {
            const Circuit a = all_kinds_circuit();
            Circuit b(3);
            for (std::size_t g = 0; g < a.size(); ++g) {
                b.add(a.gates[g]);
                if (g == 3) b.add(make_gate(K::Y, 2)).add(make_gate(K::Y, 2));
                if (g == 9) b.add(make_gate(K::S, 0)).add(make_gate(K::Sdg, 0));
            }
            CHECK(checker.check_choi(a, b).equivalent);
        }
        // Global phase uses check()'s convention A ~= global_phase * B.
        {
            const double t = 0.9;
            Circuit a(1), b(1);
            a.add(make_gate(K::RZ, 0, {t, 0, 0}));
            b.add(make_gate(K::P, 0, {t, 0, 0}));
            const EquivResult r = checker.check_choi(a, b);
            CHECK(r.equivalent);
            CHECK(std::abs(r.global_phase - std::polar(1.0, -t / 2)) < 1e-12);
            CHECK(std::abs(r.global_phase - checker.check(a, b, mixed(1, 5)).global_phase) < 1e-12);

            Circuit c = sample_circuit();
            c.global_phase = 0.7;
            const EquivResult r2 = checker.check_choi(c, sample_circuit());
            CHECK(r2.equivalent);
            CHECK(std::abs(r2.global_phase - std::polar(1.0, 0.7)) < 1e-12);
        }
    }

    void test_check_choi_not_equivalent() {
        const EquivalenceChecker checker;
        // CX vs CZ
        {
            Circuit a(2), b(2);
            a.add(make_gate(K::CX, 0, 1));
            b.add(make_gate(K::CZ, 0, 1));
            const EquivResult r = checker.check_choi(a, b);
            CHECK(!r.equivalent);
            CHECK(r.mismatch == Mismatch::State);
            CHECK(r.failing_index == 0);
            CHECK(r.failing_kind == StimulusKind::Choi);
            CHECK(r.stimuli_checked == 1);
            CHECK(r.deviation > checker.tolerance());
        }
        // Z vs I: needs a common-phase test with basis stimuli, but one Choi run catches it.
        // tr(Z I^dagger) = 0, so the output is orthogonal to |Phi>: residual sqrt(2).
        {
            Circuit z(1), id(1);
            z.add(make_gate(K::Z, 0));
            const EquivResult r = checker.check_choi(z, id);
            CHECK(!r.equivalent);
            CHECK(std::abs(r.deviation - std::sqrt(2.0)) < 1e-12);
        }
        // A tiny angle error is still caught.
        {
            Circuit a = sample_circuit(), b = sample_circuit();
            b.gates[2] = make_gate(K::RY, 2, {0.37 + 1e-6, 0, 0});
            CHECK(!checker.check_choi(a, b).equivalent);
        }
        // Conjugation matters: S vs Sdg differ, and so do the generic-matrix kinds.
        {
            Circuit a(1), b(1);
            a.add(make_gate(K::S, 0));
            b.add(make_gate(K::Sdg, 0));
            CHECK(!checker.check_choi(a, b).equivalent);
        }
        // Agrees with check() on a non-equivalent pair built from every gate kind.
        {
            const Circuit a = all_kinds_circuit();
            Circuit b = a;
            b.gates[3] = make_gate(K::X, 0);   // Y -> X
            CHECK(!checker.check_choi(a, b).equivalent);
            CHECK(!checker.check(a, b, mixed(3, 11)).equivalent);
        }
    }

    void test_check_choi_invalid_inputs() {
        const EquivalenceChecker checker;
        CHECK(throws_invalid([&] { (void)checker.check_choi(Circuit(2), Circuit(3)); }));
        CHECK(throws_invalid([&] { (void)checker.check_choi(Circuit(17), Circuit(17)); }));
        // n = 0: both circuits are scalars e^{i phase}.
        Circuit a(0), b(0);
        a.global_phase = 0.25;
        const EquivResult r = checker.check_choi(a, b);
        CHECK(r.equivalent);
        CHECK(std::abs(r.global_phase - std::polar(1.0, 0.25)) < 1e-12);
    }
} // namespace

int main() {
    test_execute_each();
    test_miter();
    test_equivalent_identities();
    test_global_phase();
    test_not_equivalent();
    test_phase_mismatch_and_lowest_index();
    test_invalid_inputs();
    test_conjugate();
    test_choi_state();
    test_choi_circuit();
    test_check_choi_equivalent();
    test_check_choi_not_equivalent();
    test_check_choi_invalid_inputs();
    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("all equiv tests passed\n");
    return 0;
}
