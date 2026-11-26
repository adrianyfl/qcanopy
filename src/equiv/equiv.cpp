#include "qcanopy/equiv.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace qcanopy {

    namespace {
        // Throws unless a Choi state over n-qubit circuits (2n qubits) fits the simulator.
        void check_choi_width(int n, const char* what) {
            if (n < 0 || 2 * n > kMaxQubits)
                throw std::invalid_argument(std::string(what) + ": 2n qubits out of range");
        }

        // The same gate moved up by `offset` qubits.
        Gate shifted(Gate g, int offset) {
            for (int i = 0; i < g.nq; ++i) g.q[i] = static_cast<uint8_t>(g.q[i] + offset);
            return g;
        }
    } // namespace

    Gate conjugate(const Gate& g) {
        const auto& p = g.params;
        switch (g.kind) {
            // Real matrices: G^* = G
            case GateKind::I: case GateKind::H: case GateKind::X: case GateKind::Z:
            case GateKind::RY:
            case GateKind::CX: case GateKind::CZ: case GateKind::SWAP:
                return g;

            // Diagonal / SX phases flip sign: conjugate = inverse
            case GateKind::S:    return make_gate(GateKind::Sdg,  g.q[0]);
            case GateKind::Sdg:  return make_gate(GateKind::S,    g.q[0]);
            case GateKind::T:    return make_gate(GateKind::Tdg,  g.q[0]);
            case GateKind::Tdg:  return make_gate(GateKind::T,    g.q[0]);
            case GateKind::SX:   return make_gate(GateKind::SXdg, g.q[0]);
            case GateKind::SXdg: return make_gate(GateKind::SX,   g.q[0]);

            // RX(t)^* = [[c, is], [is, c]] = RX(-t); RZ and P are diagonal phases
            case GateKind::RX: case GateKind::RZ: case GateKind::P:
                return make_gate(g.kind, g.q[0], {-p[0], 0, 0});

            // U(t, p, l)^* = [[c, -e^{-il} s], [e^{-ip} s, e^{-i(p+l)} c]] = U(t, -p, -l)
            case GateKind::U:
                return make_gate(GateKind::U, g.q[0], {p[0], -p[1], -p[2]});

            // Y^* = -Y has no phase-exact named form; fall through to a generic matrix.
            case GateKind::Y:
            case GateKind::U1Q:
                return make_unitary_1q(g.q[0], {std::conj(g.m[0]), std::conj(g.m[1]),
                                                std::conj(g.m[2]), std::conj(g.m[3])});
            case GateKind::U2Q: {
                std::array<cplx, 16> m{};
                for (int i = 0; i < 16; ++i) m[i] = std::conj(g.m[i]);
                return make_unitary_2q(g.q[0], g.q[1], m);
            }

            default:
                throw std::invalid_argument("conjugate: invalid GateKind");
        }
    }

    Circuit conjugate(const Circuit& c) {
        Circuit out(c.n);
        out.gates.reserve(c.size());
        for (const Gate& g : c.gates) out.add(conjugate(g));
        out.global_phase = -c.global_phase;
        return out;
    }

    StateVector choi_state(int n) {
        check_choi_width(n, "choi_state");
        StateVector state(2 * n);   // |0...0>
        const std::size_t dim = std::size_t{1} << n;
        const double amp = 1.0 / std::sqrt(static_cast<double>(dim));
        state.amplitudes[0] = {0.0, 0.0};
        for (std::size_t i = 0; i < dim; ++i) state.amplitudes[i | (i << n)] = {amp, 0.0};
        return state;
    }

    Circuit choi_prep_circuit(int n) {
        check_choi_width(n, "choi_prep_circuit");
        Circuit c(2 * n);
        for (int k = 0; k < n; ++k) {
            const auto first = static_cast<uint8_t>(k), second = static_cast<uint8_t>(k + n);
            c.add(make_gate(GateKind::H, first)).add(make_gate(GateKind::CX, first, second));
        }
        return c;
    }

    Circuit choi_circuit(const Circuit& r, const Circuit& c) {
        if (r.n != c.n) throw std::invalid_argument("choi_circuit: circuits have different qubit counts");
        check_choi_width(r.n, "choi_circuit");
        Circuit out(2 * r.n);
        out.gates.reserve(r.size() + c.size());
        for (const Gate& g : r.gates) out.add(conjugate(g));     // R^* on qubits 0..n-1
        for (const Gate& g : c.gates) out.add(shifted(g, r.n));  // C   on qubits n..2n-1
        out.global_phase = c.global_phase - r.global_phase;      // phase(R^*) = -phase(R)
        return out;
    }

    Circuit miter(const Circuit& a, const Circuit& b) {
        if (a.n != b.n) throw std::invalid_argument("miter: circuits have different qubit counts");
        Circuit m(a.n);
        m.gates.reserve(a.size() + b.size());
        m.gates = a.gates;
        for (auto it = b.gates.rbegin(); it != b.gates.rend(); ++it) m.add(inverse(*it));
        m.global_phase = a.global_phase - b.global_phase;
        return m;
    }

    EquivalenceChecker::EquivalenceChecker(double tol) : tol_(tol) {
        if (!(tol > 0.0)) throw std::invalid_argument("EquivalenceChecker: tolerance must be positive");
    }

    EquivResult EquivalenceChecker::check(const Circuit& a, const Circuit& b, const Stimulus& stimulus) const {
        if (a.n != b.n)
            throw std::invalid_argument("EquivalenceChecker::check: circuits have different qubit counts");
        if (stimulus.num_qbits != a.n)
            throw std::invalid_argument("EquivalenceChecker::check: stimulus qubit count does not match circuits");
        if (stimulus.empty())
            throw std::invalid_argument("EquivalenceChecker::check: empty stimulus");
        return check_phase_identity(miter(a, b), stimulus);
    }

    EquivResult EquivalenceChecker::check_choi(const Circuit& a, const Circuit& b) const {
        if (a.n != b.n)
            throw std::invalid_argument("EquivalenceChecker::check_choi: circuits have different qubit counts");
        check_choi_width(a.n, "EquivalenceChecker::check_choi");

        Stimulus stimulus(2 * a.n);
        stimulus.add(choi_state(a.n), StimulusKind::Choi);

        // R = A, C = B: the output is (I (x) B A^dagger)|Phi>, so the measured phase u satisfies
        // B ~= u A. Report its conjugate to match check()'s convention A ~= global_phase * B.
        EquivResult result = check_phase_identity(choi_circuit(a, b), stimulus);
        result.global_phase = std::conj(result.global_phase);
        return result;
    }

    EquivResult EquivalenceChecker::check_phase_identity(const Circuit& m, const Stimulus& stimulus) const {
        const std::size_t count = stimulus.size();

        // Per-stimulus measurements, indexed by stimulus. Written and read under `mu`.
        std::vector<cplx> phase(count);          // u_i = o_i / |o_i|
        std::vector<double> residual(count);     // || M|psi_i> - u_i |psi_i> ||
        std::vector<char> done(count, 0);        // char, not bool: distinct bytes per slot
        std::mutex mu;

        // Online early-stop test. It must only return Stop on a genuine failure under the SAME
        // rule as the final verdict below (reference phase = stimulus 0), otherwise cancelled
        // stimuli could hide an unchecked pass. Stimulus 0 is never cancelled (only indices
        // above a stopping index are), so its phase is always available eventually; results
        // that finish before it are re-checked when it arrives.
        auto on_output = [&](std::size_t i, const StateVector& out) -> Control {
            const auto& psi = stimulus.states[i].amplitudes;
            const auto& amp = out.amplitudes;

            cplx overlap{0.0, 0.0};
            for (std::size_t k = 0; k < amp.size(); ++k) overlap += std::conj(psi[k]) * amp[k];
            const double mag = std::abs(overlap);
            const cplx u = mag > 0.0 ? overlap / mag : cplx{1.0, 0.0};

            double res_sq = 0.0;
            for (std::size_t k = 0; k < amp.size(); ++k) res_sq += std::norm(amp[k] - u * psi[k]);
            const double res = std::sqrt(res_sq);

            std::lock_guard<std::mutex> lock(mu);
            phase[i] = u;
            residual[i] = res;
            done[i] = 1;

            if (res > tol_) return Control::Stop;
            if (i == 0) {
                for (std::size_t j = 1; j < count; ++j)
                    if (done[j] && residual[j] <= tol_ && std::abs(phase[j] - u) > tol_) return Control::Stop;
                return Control::Continue;
            }
            if (done[0] && std::abs(u - phase[0]) > tol_) return Control::Stop;
            return Control::Continue;
        };

        sim_.execute_each(m, stimulus, on_output);

        // Final verdict: a deterministic scan in index order over the completed stimuli,
        // so the witness is the lowest-index failure among them.
        EquivResult result;
        result.stimuli_checked = static_cast<std::size_t>(std::count(done.begin(), done.end(), 1));
        const cplx ref = phase[0];
        double max_dev = 0.0;
        for (std::size_t i = 0; i < count; ++i) {
            if (!done[i]) continue;
            const double phase_dev = std::abs(phase[i] - ref);
            Mismatch why = Mismatch::None;
            double dev = 0.0;
            if (residual[i] > tol_) { why = Mismatch::State; dev = residual[i]; }
            else if (phase_dev > tol_) { why = Mismatch::Phase; dev = phase_dev; }

            if (why != Mismatch::None) {
                result.equivalent = false;
                result.mismatch = why;
                result.deviation = dev;
                result.failing_index = i;
                result.failing_kind = stimulus.kinds[i];
                return result;
            }
            max_dev = std::max({max_dev, residual[i], phase_dev});
        }

        // No failure means no Stop was ever returned, so every stimulus ran.
        if (result.stimuli_checked != count)
            throw std::logic_error("EquivalenceChecker: stimuli cancelled without a failure");
        result.equivalent = true;
        result.global_phase = ref;
        result.deviation = max_dev;
        return result;
    }

} // namespace qcanopy
