// Tests for stimulus generation (Haar-random states) and Simulator::execute from an initial state.
#include "qcanopy/stimulus.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>

using namespace qcanopy;

namespace {
    int failures = 0;

#define CHECK(cond) \
    do { if (!(cond)) { std::printf("%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)

    double norm_sq(const StateVector& s) {
        double t = 0.0;
        for (const auto& a : s.amplitudes) t += std::norm(a);
        return t;
    }

    cplx inner(const StateVector& a, const StateVector& b) {
        cplx t{0.0, 0.0};
        for (std::size_t i = 0; i < a.amplitudes.size(); ++i) t += std::conj(a.amplitudes[i]) * b.amplitudes[i];
        return t;
    }

    bool same_state(const StateVector& a, const StateVector& b, double tol = 0.0) {
        if (a.num_qbits != b.num_qbits || a.amplitudes.size() != b.amplitudes.size()) return false;
        for (std::size_t i = 0; i < a.amplitudes.size(); ++i)
            if (std::abs(a.amplitudes[i] - b.amplitudes[i]) > tol) return false;
        return true;
    }

    void test_shape_and_norm() {
        for (int n = 0; n <= 8; ++n) {
            const StateVector s = haar_random_state(n, 42 + n);
            CHECK(s.num_qbits == n);
            CHECK(s.amplitudes.size() == (std::size_t{1} << n));
            CHECK(std::abs(norm_sq(s) - 1.0) < 1e-12);
        }
    }

    void test_determinism() {
        CHECK(same_state(haar_random_state(5, 7), haar_random_state(5, 7)));
        CHECK(!same_state(haar_random_state(5, 7), haar_random_state(5, 8), 1e-6));

        // Batch entries depend only on (n, seed, k): a longer batch extends a shorter one.
        const auto short_batch = haar_random_states(4, 3, 99);
        const auto long_batch = haar_random_states(4, 10, 99);
        CHECK(short_batch.size() == 3 && long_batch.size() == 10);
        for (std::size_t k = 0; k < short_batch.size(); ++k) CHECK(same_state(short_batch[k], long_batch[k]));
        // Batch entries are distinct from one another.
        CHECK(!same_state(long_batch[0], long_batch[1], 1e-6));
        CHECK(haar_random_states(3, 0, 1).empty());
    }

    void test_invalid_qubit_count() {
        bool threw = false;
        try { (void)haar_random_state(-1, 0); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }

    // Moments of the Haar measure on C^d: for each basis index i,
    //   E|a_i|^2 = 1/d,   E|a_i|^4 = 2/(d(d+1)),
    // and for two independent Haar states, E|<psi|phi>|^2 = 1/d.
    // Seeds are fixed, so this is deterministic; tolerances are ~5 standard errors.
    void test_haar_moments() {
        const int n = 2;
        const std::size_t d = std::size_t{1} << n;
        const std::size_t samples = 20000;
        const auto states = haar_random_states(n, samples, 2026);

        for (std::size_t i = 0; i < d; ++i) {
            double m2 = 0.0, m4 = 0.0;
            for (const auto& s : states) {
                const double p = std::norm(s.amplitudes[i]);
                m2 += p;
                m4 += p * p;
            }
            m2 /= samples;
            m4 /= samples;
            CHECK(std::abs(m2 - 1.0 / d) < 0.007);
            CHECK(std::abs(m4 - 2.0 / (d * (d + 1))) < 0.006);
        }

        // Real/imag parts should be unbiased (no preferred phase).
        cplx mean_amp{0.0, 0.0};
        for (const auto& s : states) mean_amp += s.amplitudes[0];
        mean_amp /= static_cast<double>(samples);
        CHECK(std::abs(mean_amp) < 0.02);

        double overlap = 0.0;
        for (std::size_t k = 0; k + 1 < samples; k += 2) overlap += std::norm(inner(states[k], states[k + 1]));
        overlap /= (samples / 2);
        CHECK(std::abs(overlap - 1.0 / d) < 0.01);
    }

    // Single-qubit reduced density matrix of qubit k: returns {rho00, rho01}; rho11 = 1 - rho00.
    struct Marginal { double p0; cplx coh; };
    Marginal marginal(const StateVector& s, int k) {
        const std::size_t bit = std::size_t{1} << k;
        Marginal m{0.0, {0.0, 0.0}};
        for (std::size_t i = 0; i < s.amplitudes.size(); ++i) {
            if (i & bit) continue;
            m.p0 += std::norm(s.amplitudes[i]);
            m.coh += s.amplitudes[i] * std::conj(s.amplitudes[i | bit]);
        }
        return m;
    }

    // Tr(rho_k^2): 1 iff qubit k is unentangled with the rest.
    double purity(const StateVector& s, int k) {
        const Marginal m = marginal(s, k);
        const double p1 = 1.0 - m.p0;
        return m.p0 * m.p0 + p1 * p1 + 2.0 * std::norm(m.coh);
    }

    void test_product_shape_and_determinism() {
        for (int n = 0; n <= 8; ++n) {
            const StateVector s = random_product_state(n, 300 + n);
            CHECK(s.num_qbits == n);
            CHECK(s.amplitudes.size() == (std::size_t{1} << n));
            CHECK(std::abs(norm_sq(s) - 1.0) < 1e-12);
        }
        CHECK(same_state(random_product_state(4, 5), random_product_state(4, 5)));
        CHECK(!same_state(random_product_state(4, 5), random_product_state(4, 6), 1e-6));

        const auto short_batch = random_product_states(4, 3, 99);
        const auto long_batch = random_product_states(4, 10, 99);
        for (std::size_t k = 0; k < short_batch.size(); ++k) CHECK(same_state(short_batch[k], long_batch[k]));
        CHECK(!same_state(long_batch[0], long_batch[1], 1e-6));

        // Same seed, different stimulus kind -> unrelated batches.
        CHECK(!same_state(random_product_states(1, 1, 99)[0], haar_random_states(1, 1, 99)[0], 1e-6));

        bool threw = false;
        try { (void)random_product_state(kMaxQubits + 1, 0); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }

    void test_product_is_unentangled() {
        const int n = 4;
        for (uint64_t seed = 0; seed < 20; ++seed) {
            const StateVector s = random_product_state(n, seed);
            for (int k = 0; k < n; ++k) CHECK(std::abs(purity(s, k) - 1.0) < 1e-12);
        }
        // Contrast: a Haar state is (almost surely) entangled, so the purity check has teeth.
        const StateVector h = haar_random_state(n, 0);
        for (int k = 0; k < n; ++k) CHECK(purity(h, k) < 0.99);
    }

    // Each qubit is uniform on the Bloch sphere, so p = |<0|v_k>|^2 ~ Uniform(0, 1):
    //   E[p] = 1/2, E[p^2] = 1/3,  and distinct qubits are independent: E[p_j p_k] = 1/4.
    void test_product_moments() {
        const int n = 3;
        const std::size_t samples = 20000;
        const auto states = random_product_states(n, samples, 2026);

        double m1[n] = {}, m2[n] = {}, cross01 = 0.0, cross12 = 0.0;
        for (const auto& s : states) {
            double p[n];
            for (int k = 0; k < n; ++k) {
                p[k] = marginal(s, k).p0;
                m1[k] += p[k];
                m2[k] += p[k] * p[k];
            }
            cross01 += p[0] * p[1];
            cross12 += p[1] * p[2];
        }
        for (int k = 0; k < n; ++k) {
            CHECK(std::abs(m1[k] / samples - 0.5) < 0.01);
            CHECK(std::abs(m2[k] / samples - 1.0 / 3.0) < 0.01);
        }
        CHECK(std::abs(cross01 / samples - 0.25) < 0.01);
        CHECK(std::abs(cross12 / samples - 0.25) < 0.01);
    }

    void test_execute_from_initial_state() {
        Simulator sim;
        Circuit c(3);
        c.add(make_gate(GateKind::H, 0))
         .add(make_gate(GateKind::CX, 0, 1))
         .add(make_gate(GateKind::RY, 2, {0.37, 0, 0}))
         .add(make_gate(GateKind::T, 1))
         .add(make_gate(GateKind::SWAP, 1, 2));
        c.global_phase = 0.25;

        // Starting from |000> explicitly matches the default overload.
        CHECK(same_state(sim.execute(c), sim.execute(c, StateVector(3)), 1e-12));

        // C^-1 C |psi> = |psi> for a Haar-random |psi>.
        Circuit round_trip = c;
        for (auto it = c.gates.rbegin(); it != c.gates.rend(); ++it) round_trip.add(inverse(*it));
        round_trip.global_phase = 0.0;   // phases cancel
        const StateVector psi = haar_random_state(3, 11);
        CHECK(same_state(sim.execute(round_trip, psi), psi, 1e-12));

        // Unitarity: output stays normalized and inner products are preserved.
        const StateVector phi = haar_random_state(3, 12);
        const StateVector out_psi = sim.execute(c, psi), out_phi = sim.execute(c, phi);
        CHECK(std::abs(norm_sq(out_psi) - 1.0) < 1e-12);
        CHECK(std::abs(inner(out_psi, out_phi) - inner(psi, phi)) < 1e-12);

        // Size mismatch is rejected.
        bool threw = false;
        try { (void)sim.execute(c, haar_random_state(2, 0)); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }
} // namespace

int main() {
    test_shape_and_norm();
    test_determinism();
    test_invalid_qubit_count();
    test_haar_moments();
    test_product_shape_and_determinism();
    test_product_is_unentangled();
    test_product_moments();
    test_execute_from_initial_state();
    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("all stimulus tests passed\n");
    return 0;
}
