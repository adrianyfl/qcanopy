#include "qcanopy/ir.hpp"
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace qcanopy {

    namespace {
        //          name    nq np diag   perm   selfinv sym
        constexpr GateInfo kInfo[] = {
        /* I    */ {"id",   1, 0, true,  true,  true,  false},
        /* H    */ {"h",    1, 0, false, false, true,  false},
        /* X    */ {"x",    1, 0, false, true,  true,  false},
        /* Y    */ {"y",    1, 0, false, false, true,  false},
        /* Z    */ {"z",    1, 0, true,  false, true,  false},
        /* S    */ {"s",    1, 0, true,  false, false, false},
        /* Sdg  */ {"sdg",  1, 0, true,  false, false, false},
        /* T    */ {"t",    1, 0, true,  false, false, false},
        /* Tdg  */ {"tdg",  1, 0, true,  false, false, false},
        /* SX   */ {"sx",   1, 0, false, false, false, false},
        /* SXdg */ {"sxdg", 1, 0, false, false, false, false},
        /* RX   */ {"rx",   1, 1, false, false, false, false},
        /* RY   */ {"ry",   1, 1, false, false, false, false},
        /* RZ   */ {"rz",   1, 1, true,  false, false, false},
        /* P    */ {"p",    1, 1, true,  false, false, false},
        /* U    */ {"u",    1, 3, false, false, false, false},
        /* CX   */ {"cx",   2, 0, false, true,  true,  false},
        /* CZ   */ {"cz",   2, 0, true,  false, true,  true },
        /* SWAP */ {"swap", 2, 0, false, true,  true,  true },
        /* U1Q  */ {"",     1, 0, false, false, false, false},
        /* U2Q  */ {"",     2, 0, false, false, false, false},
        };

        static_assert(sizeof(kInfo) / sizeof(kInfo[0]) == static_cast<size_t>(GateKind::Count_),
                    "kInfo must have one row per GateKind");

        constexpr double kPi = 3.14159265358979323846;
        const cplx kI{0.0, 1.0};

        cplx eix(double x) {return std::polar(1.0, x);} // e^{i x}

        void set2(Gate& g, cplx a, cplx b, cplx c, cplx d) {
            g.m = {};
            g.m[0] = a;
            g.m[1] = b;
            g.m[2] = c;
            g.m[3] = d;
        }

        constexpr double kUnitaryTol = 1e-9;

        // True if the d x d row-major matrix m satisfies M * M^dagger = I within tol.
        bool is_unitary(const cplx* m, int d, double tol) {
            for (int r = 0; r < d; ++r)
                for (int c = 0; c < d; ++c) {
                cplx sum = 0;
                for (int k = 0; k < d; ++k) sum += m[r * d + k] * std::conj(m[c * d + k]);
                const cplx expect = (r == c) ? 1.0 : 0.0;
                if (std::abs(sum - expect) > tol) return false;
                }
            return true;
        }

        int dim(const Gate& g) { return g.nq == 1 ? 2 : 4; }

        // Conjugate transpose of g's matrix, written into out (unused slots zeroed).
        std::array<cplx, 16> dagger(const Gate& g) {
            const int d = dim(g);
            std::array<cplx, 16> out{};
            for (int r = 0; r < d; ++r)
                for (int c = 0; c < d; ++c) out[r * d + c] = std::conj(g.m[c * d + r]);
            return out;
        }

        std::string gate_err(std::size_t i, const char* what) {
            return "validate: gate " + std::to_string(i) + ": " + what;
        }
    }

    const GateInfo& info(GateKind k) {
        return kInfo[static_cast<size_t>(k)];
    }

    Gate make_gate(GateKind k, uint8_t q0, std::array<double, 3> params = {0, 0, 0}) {
        Gate g = {k, 1, {q0, 0}, params};
        
        const double th = params[0];
        const double c = std::cos(th / 2), s = std::sin(th / 2);
        const double r = 1.0 / std::sqrt(2.0);

        switch (k) {
            case GateKind::I:    set2(g, 1, 0, 0, 1); break;
            case GateKind::H:    set2(g, r, r, r, -r); break;
            case GateKind::X:    set2(g, 0, 1, 1, 0); break;
            case GateKind::Y:    set2(g, 0, -kI, kI, 0); break;
            case GateKind::Z:    set2(g, 1, 0, 0, -1); break;
            case GateKind::S:    set2(g, 1, 0, 0, kI); break;
            case GateKind::Sdg:  set2(g, 1, 0, 0, -kI); break;
            case GateKind::T:    set2(g, 1, 0, 0, eix(kPi / 4)); break;
            case GateKind::Tdg:  set2(g, 1, 0, 0, eix(-kPi / 4)); break;
            case GateKind::SX:   set2(g, cplx(.5, .5), cplx(.5, -.5), cplx(.5, -.5), cplx(.5, .5)); break;
            case GateKind::SXdg: set2(g, cplx(.5, -.5), cplx(.5, .5), cplx(.5, .5), cplx(.5, -.5)); break;
            case GateKind::RX:   set2(g, c, -kI * s, -kI * s, c); break;
            case GateKind::RY:   set2(g, c, -s, s, c); break;
            case GateKind::RZ:   set2(g, eix(-th / 2), 0, 0, eix(th / 2)); break;
            case GateKind::P:    set2(g, 1, 0, 0, eix(th)); break;
            case GateKind::U: {                       // Qiskit's phase-free U(theta, phi, lambda)
                const double phi = params[1], lam = params[2];
                set2(g, c, -eix(lam) * s, eix(phi) * s, eix(phi + lam) * c);
                break;
            }
            default:
                throw std::invalid_argument("make_gate: not a named 1-qubit gate");
        }

        return g;
    }

    // for CX, CZ, SWAP
    Gate make_gate(GateKind k, uint8_t q0, uint8_t q1) {
        if (q0 == q1)
            throw std::invalid_argument("make_gate: 2-qubit gate needs two distinct qubits");

        Gate g;
        g.kind = k; g.nq = 2; g.q = {q0, q1}; g.m = {};
        auto set = [&](int row, int col, cplx v) { g.m[row * 4 + col] = v; };

        switch (k) {
            case GateKind::CX:   set(0,0,1); set(1,3,1); set(2,2,1); set(3,1,1); break;
            case GateKind::CZ:   set(0,0,1); set(1,1,1); set(2,2,1); set(3,3,-1); break;
            case GateKind::SWAP: set(0,0,1); set(1,2,1); set(2,1,1); set(3,3,1); break;
            default:
            throw std::invalid_argument("make_gate: not a named 2-qubit gate");
        }
        return g;
    }

    // U1Q
    Gate make_unitary_1q(uint8_t q0, const std::array<cplx, 4>& m) {
        if (!is_unitary(m.data(), 2, kUnitaryTol))
            throw std::invalid_argument("make_unitary_1q: matrix is not unitary");

        Gate g;
        g.kind = GateKind::U1Q; g.nq = 1; g.q = {q0, 0}; g.params = {0, 0, 0};
        g.m = {};                                  // zero all 16 slots, then copy 4
        for (int i = 0; i < 4; ++i) g.m[i] = m[i];
        return g;
    }

    // U2Q
    Gate make_unitary_2q(uint8_t q0, uint8_t q1, const std::array<cplx, 16>& m) {
        if (q0 == q1)
            throw std::invalid_argument("make_unitary_2q: needs two distinct qubits");
        if (!is_unitary(m.data(), 4, kUnitaryTol))
            throw std::invalid_argument("make_unitary_2q: matrix is not unitary");

        Gate g;
        g.kind = GateKind::U2Q; g.nq = 2; g.q = {q0, q1}; g.params = {0, 0, 0};
        g.m = m;
        return g;
    }

    Gate inverse(const Gate& g) {
        const auto& p = g.params;
        switch (g.kind) {
            // Self-inverse named gates
            case GateKind::I: case GateKind::H: case GateKind::X:
            case GateKind::Y: case GateKind::Z:
            case GateKind::CX: case GateKind::CZ: case GateKind::SWAP:
            return g;
        
            // Named inverse pairs
            case GateKind::S:    return make_gate(GateKind::Sdg,  g.q[0]);
            case GateKind::Sdg:  return make_gate(GateKind::S,    g.q[0]);
            case GateKind::T:    return make_gate(GateKind::Tdg,  g.q[0]);
            case GateKind::Tdg:  return make_gate(GateKind::T,    g.q[0]);
            case GateKind::SX:   return make_gate(GateKind::SXdg, g.q[0]);
            case GateKind::SXdg: return make_gate(GateKind::SX,   g.q[0]);
        
            // Rotations: negate the angle
            case GateKind::RX: case GateKind::RY: case GateKind::RZ: case GateKind::P:
                return make_gate(g.kind, g.q[0], {-p[0], 0, 0});
        
            // U(theta, phi, lambda)^dagger = U(-theta, -lambda, -phi)
            case GateKind::U:
                return make_gate(GateKind::U, g.q[0], {-p[0], -p[2], -p[1]});
        
            // Generic matrices: conjugate transpose
            case GateKind::U1Q: {
                const auto d = dagger(g);
                return make_unitary_1q(g.q[0], {d[0], d[1], d[2], d[3]});
            }
            case GateKind::U2Q:
                return make_unitary_2q(g.q[0], g.q[1], dagger(g));
        
            default:
                throw std::invalid_argument("inverse: invalid GateKind");
        }
    }

    void validate(const Circuit& c) {
        if (c.n < 1 || c.n > kMaxQubits)
            throw std::invalid_argument("validate: qubit count " + std::to_string(c.n) +
                                        " outside [1, " + std::to_string(kMaxQubits) + "]");
        if (!std::isfinite(c.global_phase))
            throw std::invalid_argument("validate: global_phase is not finite");
        
        for (std::size_t i = 0; i < c.gates.size(); ++i) {
            const Gate& g = c.gates[i];
        
            if (static_cast<size_t>(g.kind) >= static_cast<size_t>(GateKind::Count_))
                throw std::invalid_argument(gate_err(i, "invalid kind"));
            if (g.nq != 1 && g.nq != 2)
                throw std::invalid_argument(gate_err(i, "nq must be 1 or 2"));
            if (g.nq != info(g.kind).nq)
                throw std::invalid_argument(gate_err(i, "nq does not match gate kind"));
        
            for (int j = 0; j < g.nq; ++j)
                if (g.q[j] >= c.n)
                    throw std::invalid_argument(gate_err(i, "qubit index out of range"));
            if (g.nq == 2 && g.q[0] == g.q[1])
                throw std::invalid_argument(gate_err(i, "2-qubit gate acts on the same qubit twice"));
        
            const int d = dim(g);
            for (int j = 0; j < d * d; ++j)
                if (!std::isfinite(g.m[j].real()) || !std::isfinite(g.m[j].imag()))
                    throw std::invalid_argument(gate_err(i, "matrix has non-finite entries"));
        }
    }

    bool same_gate(const Gate& a, const Gate& b, double tol) {
        if (a.kind != b.kind || a.nq != b.nq) return false;
        
        if (a.nq == 1) {
            if (a.q[0] != b.q[0]) return false;
        } else {
            const bool same_order = a.q[0] == b.q[0] && a.q[1] == b.q[1];
            const bool swapped    = a.q[0] == b.q[1] && a.q[1] == b.q[0];
            // Symmetric gates (CZ, SWAP) have the same matrix in either qubit order.
            if (!(same_order || (swapped && info(a.kind).symmetric))) return false;
        }
        
        const int d = dim(a);
        for (int j = 0; j < d * d; ++j)
            if (std::abs(a.m[j] - b.m[j]) > tol) return false;
        return true;
    }
}
