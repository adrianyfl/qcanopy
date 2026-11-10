// qcanopy/peel.hpp — strip the common prefix and suffix off two circuits.
//
// Given a reference circuit R and a candidate C, peel() removes the gates both circuits
// start with and end with, leaving only the parts that differ (R_mid, C_mid). The peeled
// gates are not kept: the equivalence check only has to look at the miter C_mid · R_mid†
// (see 3), so the rewrite loop calls peel() until it removes nothing.
//
// 1. Frontier peeling. A gate is peelable when it is at the front (or back) of every wire it
//    touches, in both circuits. On each wire the front gate of R is compared with the front
//    gate of C; if they match, both are removed, which may expose new front gates. Each wire
//    has one front gate per circuit, so every match is forced and the result is the unique
//    maximal common prefix (and suffix). Gates are never commuted; that is the rewriter's job.
//
// 2. Matching up to phase. Two gates match if they act on the same qubits (same order) and
//    their matrices satisfy C_gate = e^{i alpha} R_gate. Gate kinds are not compared, so a
//    U1Q can match a U or H with the same matrix.
//
// 3. Phase bookkeeping. Each alpha is added to C.global_phase. A global phase commutes with
//    everything, so C · R† = S · (C_mid · R_mid†) · S† with the phases tracked exactly:
//    the miter is the identity (up to phase) iff R and C are equivalent, and
//    C.global_phase - R.global_phase is the phase by which they differ.
//
// Run normalize() on both DAGs first so that equal gates are spelled the same way.

#pragma once

#include "qcanopy/dag.hpp"

namespace qcanopy {

    // Removes the maximal common prefix, then the maximal common suffix, from R and C.
    // Returns the number of gate pairs removed (0 means nothing matched). R and C must have
    // the same number of qubits; throws std::invalid_argument otherwise. On return they hold
    // R_mid and C_mid; removed nodes are tombstoned (alive = false), not compacted.
    std::size_t peel(CircuitDAG& R, CircuitDAG& C);

}
