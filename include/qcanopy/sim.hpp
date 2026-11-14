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

    // Which generator produced a stimulus state (see qcanopy/stimulus.hpp).
    enum class StimulusKind : uint8_t {
        Custom,     // added by hand
        Haar,       // Haar-random state
        Product,    // random product state
        Basis,      // computational basis state
    };

    // A batch of input states for one circuit width, possibly mixing several kinds.
    // Invariant: every state has num_qbits qubits, and kinds[i] describes states[i].
    // Each state is simulated independently, so a batch is the unit of parallel work.
    struct Stimulus {
        int num_qbits;
        std::vector<StateVector> states;
        std::vector<StimulusKind> kinds;

        explicit Stimulus(int n) : num_qbits(n) {}

        // Appends one state. Throws std::invalid_argument if s.num_qbits != num_qbits.
        Stimulus& add(StateVector s, StimulusKind kind = StimulusKind::Custom);

        // Appends all of `other`, preserving order. Throws std::invalid_argument on width mismatch.
        Stimulus& append(const Stimulus& other);

        std::size_t size() const { return states.size(); }
        bool empty() const { return states.empty(); }
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

        // Executes the circuit on every stimulus state; result[i] is the output for
        // stimulus.states[i]. Throws std::invalid_argument if stimulus.num_qbits != circuit.n.
        std::vector<StateVector> execute(const Circuit& circuit, const Stimulus& stimulus) const;

        // Extracts the topological schedule from the DAG and executes it.
        StateVector execute(const CircuitDAG& dag) const;

    private:
        // Core bitwise application of unitary matrices.
        void apply_1q_gate(StateVector& state, const Gate& gate) const;
        void apply_2q_gate(StateVector& state, const Gate& gate) const;
    };

} // namespace qcanopy