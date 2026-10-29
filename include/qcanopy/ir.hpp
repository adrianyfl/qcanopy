// qcanopy/ir.hpp — the circuit intermediate representation (IR).
//
// 1. Qubit ordering: qubit 0 is the LEAST significant bit of an amplitude index
//    (same as Qiskit). Amplitude index i = b_{n-1} ... b_1 b_0, b_k = qubit k.
//
// 2. Gate matrices are stored ROW-MAJOR in Gate::m, always in double precision.
//      1-qubit gate: m[0..3]  = [[m0, m1], [m2, m3]]
//      2-qubit gate: m[0..15] = 4x4 matrix
//
// 3. Local basis for 2-qubit gates: the local index is (b_q[1] << 1) | b_q[0],
//    i.e. q[0] is the low bit (again Qiskit's convention for qargs [q0, q1]).
//    For controlled gates q[0] is the CONTROL and q[1] is the TARGET, so
//    CX(control=q[0], target=q[1]) = [[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]].
//
// 4. Matrices are exact unitaries, including global phase. Equivalence checks
//    ignore global phase; the IR does not.

#pragma once

#include <array>
#include <complex>
#include <cstdint>
#include <string_view>
#include <vector>

namespace qcanopy {
    
    using cplx = std::complex<double>;

    // Maximum qubits per circuit
    inline constexpr int kMaxQubits = 32;

    // Gate Kinds
    enum class GateKind : uint8_t {
        // 1-qubit fixed
        I, H, X, Y, Z, S, Sdg, T, Tdg, SX, SXdg,
        // 1-qubit parameterized
        RX, RY, RZ, P,  // P(lambnda) = diag(1, e^{i lambda})
        U,              // U(theta, phi, lambda)
        // 2-qubit fixed
        CX, CZ, SWAP,
        // Generic unitaries produced by rewriter, Matrix only.
        U1Q, U2Q,
        Count_          // number of kinds
    };

    // Gate Information
    struct GateInfo {
        std::string_view name;  // lower-case OpenQASM name ("" for U1Q/U2Q)
        uint8_t nq;             // number of qubits
        uint8_t nparams;        // number of angle parameters 
        bool diagonal;          // matrix is diagonal (only phases)
        bool permutation;       // matrix is a 0/1 permutation (only swap amplitudes)
        bool self_inverse;      // G * G = I
        bool symmetric;         // 2-qubit gate unchanged by swapping q[0], q[1]
    };

    const GateInfo& info(GateKind k);   // table lookup

    // Gate
    struct Gate {
        GateKind kind = GateKind::I;
        uint8_t nq = 1;
        std::array<uint8_t, 2> q{0, 0};         // qubits
        std::array<double, 3> params{0, 0, 0};  // angles
        std::array<cplx, 16> m{};               // matrix, row-major

        const cplx& at(int r, int c) const {return m[r * (nq == 1 ? 2 : 4) + c];}
        cplx& at(int r, int c) {return m[r * (nq == 1 ? 2 : 4) + c];}
    };

    // Constructors for gates
    Gate make_gate(GateKind k, uint8_t q0, std::array<double, 3> params = {0, 0, 0});
    Gate make_gate(GateKind k, uint8_t q0, uint8_t q1);                                 // for CX, CZ, SWAP
    Gate make_unitary_1q(uint8_t q0, const std::array<cplx, 4>& m);                     // U1Q
    Gate make_unitary_2q(uint8_t q0, uint8_t q1, const std::array<cplx, 16>& m);        // U2Q

    // Inverse gate, return U1Q/U2Q if no native inverse
    Gate inverse(const Gate& g);

    // True if the two gates act on any common qubit
    inline bool overlaps(const Gate& a, const Gate& b) {
        for (int i = 0; i < a.nq; ++i)
            for (int j = 0; j < b.nq; ++j)
                if (a.q[i] == b.q[j]) return true;
        return false;
    }

    // Circuit
    struct Circuit {
        int n = 0;                  // number of qubits
        std::vector<Gate> gates;    // applied in order
        double global_phase = 0.0;    // to support OpenQASM 3.x

        Circuit() = default;
        explicit Circuit(int nqubits) : n(nqubits) {}

        Circuit& add(const Gate& g) {gates.push_back(g); return *this;}
        std::size_t size() const {return gates.size();}
    };

    // Checks if a circuit is valid, throws std::invalid_argument if the circuit is malformed
    void validate(const Circuit& c);

    // Gate-level structural equality
    bool same_gate(const Gate& a, const Gate& b, double tol = 1e-12);
}
