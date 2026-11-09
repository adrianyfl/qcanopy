#include "qcanopy/sim.hpp"

namespace qcanopy {
    // state-vector representation
    class StateVector {
    public:
        int n_qubits;
        std::vector<std::complex<double>> amplitudes;

        // Initializes state vector.
        explicit StateVector(int n) : n_qubits(n){
            // caculate state size
            size_t total_states = 1ULL << n;
            // sets every amplitude to 0.0 + 0.0i
            amplitudes.assign(total_states, {0.0, 0.0});
            // We give |0...0> an amplitude of 1.0 + 0.0i
            amplitudes[0] = {1.0, 0.0};
        }

        // Returns the probability of a state
        double probability(size_t state_idx) const{
            return std::norm(amplitudes[state_idx]);
        }

        // Applies a global phase to the entire state vector.
        void apply_global_phase(double phase){
            // If the phase is 0, skip 
            if (std::abs(phase) < 1e-9) {
                return;
            }
            // Calculate the complex phase factor: e^{i * phase}
            std::complex<double> phase_factor = std::polar(1.0, phase);
            // multiply each amplitude with phase factor.
            for (auto& amp : amplitudes) {
                amp[i] *= phase_factor;
            }
        }
    };
    // Core simulator engine.
    class Simulator {
    public:
        Simulator() = default;

        // Executes a flattened circuit and returns the final state.
        Simulator::execute(const Circuit& circuit) {
            // Initialize the state vector
            StateVector state(circuit.n);

            // Iterate through the gates and apply them
            for (const auto& gate : circuit.gates) {
                if (gate.nq == 1) {
                    apply_1q_gate(state, gate);
                } else if (gate.nq == 2) {
                    apply_2q_gate(state, gate);
                } else {
                    throw std::runtime_error("Simulator only supports 1-qubit and 2-qubit gates.");
                }
            }
            // Return final state vector
            return state;
        }

        // Execute directly from the DAG
        StateVector Simulator::execute(const CircuitDAG& dag) {
            // Flatten the graph into a linear schedule, then pass it to the main execute function
            StateVector final_state = execute(dag.to_circuit());
            
            return final_state;
        }

    private:
        
        // Core bitwise application of unitary matrices.
        void apply_1q_gate(StateVector& state, const Gate& gate) const{
            return
        }
        void apply_2q_gate(StateVector& state, const Gate& gate) const{
            return
        }
    };
} // namespace qcanopy