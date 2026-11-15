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

        double tolerance() const { return tol_; }
        const Simulator& simulator() const { return sim_; }

    private:
        Simulator sim_;
        double tol_;
    };

    // The miter B^-1 A: A's gates, then B's gates inverted in reverse order, with global
    // phase a.global_phase - b.global_phase. Equals e^{i phi} I iff A = e^{i phi} B.
    // Throws std::invalid_argument if a.n != b.n.
    Circuit miter(const Circuit& a, const Circuit& b);

} // namespace qcanopy
