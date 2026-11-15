#pragma once

#include <atomic>
#include <cstddef>
#include <functional>
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

    // Returned by an execute_each() callback: keep going, or cancel the remaining states.
    enum class Control : uint8_t { Continue, Stop };

    // Callback for execute_each(): receives (stimulus index, output state).
    using OutputVisitor = std::function<Control(std::size_t, const StateVector&)>;

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

        // Streams the batch: runs the circuit on each stimulus state and calls
        // on_output(i, output) as soon as state i finishes; the output is discarded afterwards.
        // If on_output returns Control::Stop for index i, every state with index > i is
        // cancelled (skipped if not started, abandoned at the next checkpoint if running).
        // States with index < i always run to completion, so the caller sees every result
        // up to the first stop: the lowest failing index is found regardless of timing.
        //
        // Currently serial and in index order, but callers must not rely on that: once this
        // is parallelized, on_output may be called concurrently and out of order, so it must
        // be thread-safe. Throws std::invalid_argument if stimulus.num_qbits != circuit.n.
        void execute_each(const Circuit& circuit, const Stimulus& stimulus, const OutputVisitor& on_output) const;

        // Extracts the topological schedule from the DAG and executes it.
        StateVector execute(const CircuitDAG& dag) const;

    private:
        // Gates applied between cancellation checks inside a single run.
        static constexpr std::size_t kCancelCheckInterval = 64;

        // Applies the circuit (gates + global phase) to `state` in place. If stop_at is non-null,
        // abandons the run as soon as index > *stop_at, checked every kCancelCheckInterval gates.
        // Returns false if the run was abandoned (state is then partially evolved).
        bool evolve(const Circuit& circuit, StateVector& state,
                    const std::atomic<std::size_t>* stop_at = nullptr, std::size_t index = 0) const;

        // Core bitwise application of unitary matrices.
        void apply_1q_gate(StateVector& state, const Gate& gate) const;
        void apply_2q_gate(StateVector& state, const Gate& gate) const;
    };

} // namespace qcanopy