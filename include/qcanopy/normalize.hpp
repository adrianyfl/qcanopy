// qcanopy/normalize.hpp — put a circuit DAG into canonical form before rewriting.
//
// Two circuits that differ only in how they spell the same gates should normalize to the same
// DAG, so that later passes see their shared prefixes and suffixes. normalize() preserves the
// exact unitary, including global phase (tracked in CircuitDAG::global_phase), up to kAngleTol.
//
// 1. Rotation families are written one way. Each 1-qubit gate about a fixed axis is read as
//    e^{i a} R(theta) and rewritten with theta wrapped into (-pi, pi]:
//      Z axis (I, Z, S, Sdg, T, Tdg, RZ, P, U(0, phi, lambda)) -> Z/S/Sdg/T/Tdg if theta is a
//                                                                named angle, else P(theta)
//      X axis (X, SX, SXdg, RX)                               -> X/SX/SXdg, else RX(theta)
//      Y axis (Y, RY)                                         -> Y, else RY(theta)
//    Identity rotations are removed. The phase a goes into global_phase.
// 2. U(theta, phi, lambda) has all three angles wrapped into (-pi, pi] and theta >= 0.
// 3. Symmetric 2-qubit gates (CZ, SWAP) have q[0] < q[1].
// 4. Neighbouring gates on a wire from the same rotation family are fused into one.
// 5. Neighbouring inverse pairs are cancelled (2-qubit pairs must be neighbours on both wires).
//
// Steps 4 and 5 repeat until nothing changes, so the result has no fusable or cancellable
// neighbours. Gates are never moved past each other; commutation is left to the rewriter.

#pragma once

#include "qcanopy/dag.hpp"

namespace qcanopy {

    // Angles within this distance of a named angle are snapped to it (and of 0, dropped).
    inline constexpr double kAngleTol = 1e-10;

    // Normalizes dag in place. Removed nodes are tombstoned (alive = false), not compacted.
    void normalize(CircuitDAG& dag);

}
