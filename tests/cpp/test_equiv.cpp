// Tests for EquivalenceChecker, miter(), and Simulator::execute_each (streaming + early stop).
#include "qcanopy/equiv.hpp"
#include "qcanopy/stimulus.hpp"

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
} // namespace

int main() {
    test_execute_each();
    test_miter();
    test_equivalent_identities();
    test_global_phase();
    test_not_equivalent();
    test_phase_mismatch_and_lowest_index();
    test_invalid_inputs();
    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("all equiv tests passed\n");
    return 0;
}
