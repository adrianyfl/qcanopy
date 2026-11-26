// qcanopy/equiv.hpp — simulation-based equivalence checking.
//
// A and B are equivalent iff A = e^{i phi} B for one global phase phi (ir.hpp, point 4).
// The checker owns a Simulator and drives it with a Stimulus (qcanopy/stimulus.hpp):
//
//   1. Build the miter M = B^-1 A once. If A ~ B then M = e^{i phi} I.
//   2. For each stimulus |psi_i>, simulate M|psi_i> and compute the overlap
//      o_i = <psi_i|M|psi_i> and its phase u_i = o_i / |o_i|.
//   3. The pair passes iff, for every i,
//        residual_i = || M|psi_i> - u_i |psi_i> ||  <= tol   (output is a phase times the input)
//        | u_i - u_ref |                           <= tol   (that phase is the SAME for all i)
//      where u_ref is the phase of the lowest-index stimulus. The common-phase test is what
//      catches Z vs I on basis inputs (see stimulus.hpp).
//   4. The first failure stops the simulator early (Simulator::execute_each); outputs are
//      compared as they stream out and never stored.
//
// One miter run per stimulus needs one state vector instead of two (A and B separately)
// for the same total gate count.
//
// Interpreting the verdict: NotEquivalent is definitive (a witness stimulus is reported).
// Equivalent means no stimulus distinguished the circuits. With Haar-random stimuli a
// non-equivalent pair slips through with probability 0 (up to tol); with all_basis_states(n)
// it is a proof (up to tol).
//
// ---- Choi-state checking (check_choi) ----
//
// An alternative that needs no stimulus and gives a proof in ONE simulation, at the cost of
// doubling the width to 2n qubits. Take the maximally entangled (Choi) state
//
//   |Phi> = 2^{-n/2} sum_i |i>_first |i>_second
//
// where the first half is qubits 0..n-1 and the second half is qubits n..2n-1. Apply R^*
// (the complex conjugate of R, NOT the inverse) to the first half and C to the second half.
// Using the transpose trick (M (x) I)|Phi> = (I (x) M^T)|Phi> with (R^*)^T = R^dagger:
//
//   (R^* (x) C)|Phi> = (I (x) C R^dagger)|Phi>
//
// so the output is a phase times |Phi> iff C R^dagger = e^{i phi} I, i.e. iff C ~ R. The overlap
// <Phi|out> = tr(C R^dagger) / 2^n, whose squared magnitude is the process fidelity, and the
// residual || out - u|Phi> || equals the normalized Frobenius distance ||C R^dagger - u I||_F / 2^{n/2}.
// The same pass/fail test as above then applies to this single 2n-qubit stimulus.
//
// Cost: one run of |A| + |B| gates on a 4^n-amplitude state (16 * 4^n bytes: n = 12 -> 256 MiB),
// versus 2^n runs of 2^n amplitudes for all_basis_states(n). Requires 2n <= kMaxQubits.

#pragma once

#include "qcanopy/ir.hpp"
#include "qcanopy/dag.hpp"
#include "qcanopy/sim.hpp"
#include <complex>
#include <cstddef>
#include <cstdint>

namespace qcanopy {

    // Which test a failing stimulus violated.
    enum class Mismatch : uint8_t {
        None,       // no failure
        State,      // M|psi> is not a phase multiple of |psi>
        Phase,      // M|psi> = u|psi>, but u differs from the reference phase
    };

    struct EquivResult {
        bool equivalent = false;
        std::size_t stimuli_checked = 0;    // outputs compared before stopping

        // If equivalent: A ~= global_phase * B, measured on the lowest-index stimulus.
        cplx global_phase{1.0, 0.0};

        // Equivalent:     the largest deviation seen (<= tol).
        // Not equivalent: the deviation at the failing stimulus (> tol).
        double deviation = 0.0;

        // Witness (valid iff !equivalent): the lowest-index failing stimulus.
        Mismatch mismatch = Mismatch::None;
        std::size_t failing_index = SIZE_MAX;
        StimulusKind failing_kind = StimulusKind::Custom;
    };

    class EquivalenceChecker {
    public:
        // tol bounds both the residual 2-norm and the phase difference (both scale-free,
        // since states are normalized). Throws std::invalid_argument if tol <= 0.
        explicit EquivalenceChecker(double tol = 1e-9);

        // Checks A against B on every state of `stimulus`, stopping at the first failure.
        // Throws std::invalid_argument if a.n != b.n, stimulus.num_qbits != a.n,
        // or the stimulus is empty (an empty check would pass vacuously).
        EquivResult check(const Circuit& a, const Circuit& b, const Stimulus& stimulus) const;

        // Checks A against B with the Choi state: simulates choi_circuit(a, b), i.e. A^* on the
        // first half and B on the second half, on choi_state(a.n) and tests that the output is a
        // phase times the input. A single run, and a proof (up to tol) either way.
        // The result reads exactly like check(): global_phase is phi with A ~= e^{i phi} B,
        // stimuli_checked is 1, and a failure is Mismatch::State at index 0, kind Choi.
        // Throws std::invalid_argument if a.n != b.n or 2 * a.n > kMaxQubits.
        EquivResult check_choi(const Circuit& a, const Circuit& b) const;

        double tolerance() const { return tol_; }
        const Simulator& simulator() const { return sim_; }

    private:
        // Runs `m` on every stimulus state and tests m|psi_i> = u|psi_i> with one common u
        // (the algorithm described at the top of this file). Shared by check() and check_choi().
        EquivResult check_phase_identity(const Circuit& m, const Stimulus& stimulus) const;

        Simulator sim_;
        double tol_;
    };

    // The miter B^-1 A: A's gates, then B's gates inverted in reverse order, with global
    // phase a.global_phase - b.global_phase. Equals e^{i phi} I iff A = e^{i phi} B.
    // Throws std::invalid_argument if a.n != b.n.
    Circuit miter(const Circuit& a, const Circuit& b);

    // ---- Choi-state building blocks ----

    // The elementwise complex conjugate G^* of a gate (not the inverse). Named gates map to
    // named gates where one exists (S -> Sdg, T -> Tdg, SX -> SXdg, RX/RZ/P(t) -> (-t),
    // U(t, p, l) -> U(t, -p, -l), real gates unchanged); Y and generic matrices become U1Q/U2Q.
    Gate conjugate(const Gate& g);

    // The complex conjugate circuit C^*: every gate conjugated, global phase negated.
    Circuit conjugate(const Circuit& c);

    // The Choi state |Phi> = 2^{-n/2} sum_i |i>|i> on 2n qubits: qubit k (first half) is paired
    // with qubit k + n (second half), so amplitude index i | (i << n) holds 2^{-n/2}.
    // Throws std::invalid_argument if n < 0 or 2n > kMaxQubits.
    StateVector choi_state(int n);

    // A 2n-qubit circuit preparing choi_state(n) from |0...0>: H on qubit k, then CX(k, k + n),
    // for k = 0..n-1. Useful when the Choi state should come from the simulator itself.
    // Throws std::invalid_argument if n < 0 or 2n > kMaxQubits.
    Circuit choi_prep_circuit(int n);

    // The 2n-qubit circuit R^* (x) C: R's gates conjugated on qubits 0..n-1, C's gates shifted to
    // qubits n..2n-1, global phase c.global_phase - r.global_phase. Applied to choi_state(n) it
    // yields (I (x) C R^dagger)|Phi>. Throws std::invalid_argument if r.n != c.n or 2n > kMaxQubits.
    Circuit choi_circuit(const Circuit& r, const Circuit& c);

} // namespace qcanopy
