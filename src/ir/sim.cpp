#include "qcanopy/sim.hpp"

namespace qcanopy {
    // state-vector representation
    class StateVector {
    public:
        int n_qbits;
        std::vector<cplx> amplitudes;

        // Initializes state vector.
        explicit StateVector(int n) : n_qbits(n){
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
            cplx phase_factor = std::polar(1.0, phase);
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
            int q = gate.q[0];                   // The target qubit
            size_t bit_mask = 1ULL << q;         // A mask for the q-th bit
            size_t half_dim = state.amplitude.size() / 2;

            // Grab the matrix elements
            cplx u00 = gate.m[0], u01 = gate.m[1];
            cplx u10 = gate.m[2], u11 = gate.m[3];

            // Loop 2^(n-1) times
            for (size_t i = 0; i < half_dim; ++i) {
                 
                // grab the index pair s.t. only the target qbit differ
                // For example (target q-th bit) q = 1, n_qbit = 3
                // when i = 00 => idx0 = 000, idx1 = 010
                // when i = 01 => idx0 = 001, idx1 = 011
                // when i = 10 => idx0 = 100, idx1 = 110
                // when i = 11 => idx0 = 101, idx1 = 111

                // Splits the binary number 'i' in half, shifts the top half left by 1, 
                // and glues the bottom half back on.
                size_t idx0 = ((i >> q) << (q + 1)) | (i & (bit_mask - 1));
                
                // The paired index, with the q-th bit set to 1.
                size_t idx1 = idx0 | bit_mask;

                // Fetch the two complex amplitudes
                cplx a0 = state.amplitudes[idx0];
                cplx a1 = state.amplitudes[idx1];

                // Apply the 2x2 matrix multiplication directly to this pair!
                state.amplitudes[idx0] = u00 * a0 + u01 * a1;
                state.amplitudes[idx1] = u10 * a0 + u11 * a1;
            }
        }
        void apply_2q_gate(StateVector& state, const Gate& gate) const{
            int q0 = gate.q[0]; // control qubit
            int q1 = gate.q[1]; // target qubit

            // We must insert the zeros from right to left (lowest bit first) 
            // so we don't mess up the offset for the second insertion.
            int q_min = std::min(q0, q1);
            int q_max = std::max(q0, q1);

            size_t mask0 = 1ULL << q0;
            size_t mask1 = 1ULL << q1;

            size_t quarter_dim = state.amplitudes.size() / 4;
            const auto& m = gate.m; // grab gate matrix

            for (size_t i = 0; i < quarter_dim; ++i) {
                
                // Insert 0 at the q_min position
                size_t temp = ((i >> q_min) << (q_min + 1)) | (i & ((1ULL << q_min) - 1));
                
                // Insert 0 at the q_max position
                size_t idx00 = ((temp >> q_max) << (q_max + 1)) | (temp & ((1ULL << q_max) - 1));

                // Produce the other 3 pairs
                size_t idx01 = idx00 | mask1;           // q0 = 0, q1 = 1
                size_t idx10 = idx00 | mask0;           // q0 = 1, q1 = 0
                size_t idx11 = idx00 | mask0 | mask1;   // q0 = 1, q1 = 1

                // Grab amplitudes
                cplx a00 = state.amplitudes[idx00];
                cplx a01 = state.amplitudes[idx01];
                cplx a10 = state.amplitudes[idx10];
                cplx a11 = state.amplitudes[idx11];

                // Apply 4x4 matrix multiplication
                state.amplitudes[idx00] = m[0]*a00  + m[1]*a01  + m[2]*a10  + m[3]*a11;
                state.amplitudes[idx01] = m[4]*a00  + m[5]*a01  + m[6]*a10  + m[7]*a11;
                state.amplitudes[idx10] = m[8]*a00  + m[9]*a01  + m[10]*a10 + m[11]*a11;
                state.amplitudes[idx11] = m[12]*a00 + m[13]*a01 + m[14]*a10 + m[15]*a11;
            }
        }
    }
} // namespace qcanopy