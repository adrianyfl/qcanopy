// qcanopy/stimulus.hpp — input states ("stimuli") for simulation-based equivalence checking.
//
// Two circuits A and B are equivalent (up to global phase) iff A|psi> ~ B|psi> for every |psi>.
// A checker drives both circuits with a set of stimulus states and compares the outputs; this
// header produces those states.
//
// Determinism: every generator here is a pure function of (n, seed). Sampling uses an in-house
// SplitMix64 + Box-Muller pipeline rather than <random> distributions, whose output differs
// between standard libraries (MSVC vs libstdc++), so stimuli are reproducible across platforms.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "qcanopy/sim.hpp"

namespace qcanopy {

    // A Haar-random pure state on n qubits: 2^n i.i.d. standard complex Gaussians, normalized.
    // The resulting vector is uniformly distributed on the unit sphere of C^(2^n), which is
    // exactly the unitarily invariant (Haar) measure on pure states.
    // Throws std::invalid_argument if n is outside [0, kMaxQubits].
    StateVector haar_random_state(int n, uint64_t seed);

    // `count` independent Haar-random states. State k depends only on (n, seed, k), so
    // increasing `count` extends the batch without changing existing entries.
    std::vector<StateVector> haar_random_states(int n, std::size_t count, uint64_t seed);

    // A random product state |v_{n-1}> (x) ... (x) |v_1> (x) |v_0>, where each single-qubit
    // state |v_k> is independently Haar-random (uniform on the Bloch sphere). Unentangled, so it
    // models inputs a real device can prepare with one layer of 1-qubit gates.
    // Throws std::invalid_argument if n is outside [0, kMaxQubits].
    StateVector random_product_state(int n, uint64_t seed);

    // `count` independent random product states; same batch semantics as haar_random_states.
    // Batches use a separate seed domain, so the same seed yields unrelated Haar and product batches.
    std::vector<StateVector> random_product_states(int n, std::size_t count, uint64_t seed);

} // namespace qcanopy
