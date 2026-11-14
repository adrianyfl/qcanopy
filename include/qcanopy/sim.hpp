#pragma once

#include <vector>
#include <complex>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include "qcanopy/ir.hpp"
#include "qcanopy/dag.hpp"

namespace qcanopy {
    // A simple state-vector representation of a quantum system.
    class StateVector {
    public:
        int num_qbits;
        std::vector<cplx> amplitudes;

        // Initializes the |0...0> state for n qubits.
        explicit StateVector(int n);

        // Returns the probability of measuring a specific basis state (0 to 2^n - 1).
        double probability(size_t state_index) const;

        // Applies a global phase to the entire state vector.
        void apply_global_phase(double phase);
    };

    // Core simulator engine.
    class Simulator {
    public:
        Simulator() = default;

        // Executes a flattened circuit and returns the final state.
        StateVector execute(const Circuit& circuit) const;

        // Executes a flattened circuit starting from `initial` instead of |0...0>.
        // Throws std::invalid_argument if initial.num_qbits != circuit.n.
        StateVector execute(const Circuit& circuit, StateVector initial) const;

        // Extracts the topological schedule from the DAG and executes it.
        StateVector execute(const CircuitDAG& dag) const;

    private:
        // Core bitwise application of unitary matrices.
        void apply_1q_gate(StateVector& state, const Gate& gate) const;
        void apply_2q_gate(StateVector& state, const Gate& gate) const;
    };

} // namespace qcanopy