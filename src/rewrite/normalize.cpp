#include "qcanopy/normalize.hpp"

#include <cmath>
#include <vector>

namespace qcanopy {

    namespace {
        constexpr double kPi = 3.14159265358979323846;
        constexpr double kTwoPi = 2 * kPi;
        constexpr uint32_t kNone = CircuitDAG::kNone;

        bool near(double a, double b) { return std::abs(a - b) < kAngleTol; }

        // Wraps x into (-pi, pi] and returns k such that x_in = x_out + 2*pi*k.
        double wrap(double& x) {
            double k = std::round(x / kTwoPi);
            x -= kTwoPi * k;
            if (x <= -kPi + kAngleTol) { x += kTwoPi; k -= 1; }
            return k;
        }

        enum class Axis : uint8_t { None, X, Y, Z };

        // If g is a rotation about a fixed axis, sets g = e^{i phase} R(angle) and returns the axis.
        // R is RX for X, RY for Y, and P for Z.
        Axis as_rotation(const Gate& g, double& angle, double& phase) {
            phase = 0;
            switch (g.kind) {
                case GateKind::I:    angle = 0;        return Axis::Z;
                case GateKind::Z:    angle = kPi;      return Axis::Z;
                case GateKind::S:    angle = kPi / 2;  return Axis::Z;
                case GateKind::Sdg:  angle = -kPi / 2; return Axis::Z;
                case GateKind::T:    angle = kPi / 4;  return Axis::Z;
                case GateKind::Tdg:  angle = -kPi / 4; return Axis::Z;
                case GateKind::P:    angle = g.params[0]; return Axis::Z;
                case GateKind::RZ:                                        // RZ(t) = e^{-it/2} P(t)
                    angle = g.params[0]; phase = -angle / 2; return Axis::Z;
                case GateKind::U: {                                       // U(2pi k, phi, lam) = (-1)^k P(phi + lam)
                    double theta = g.params[0];
                    const double k = wrap(theta);
                    if (!near(theta, 0)) return Axis::None;
                    angle = g.params[1] + g.params[2]; phase = k * kPi; return Axis::Z;
                }
                case GateKind::X:    angle = kPi;      phase = kPi / 2;  return Axis::X;  // X = i RX(pi)
                case GateKind::SX:   angle = kPi / 2;  phase = kPi / 4;  return Axis::X;  // SX = e^{i pi/4} RX(pi/2)
                case GateKind::SXdg: angle = -kPi / 2; phase = -kPi / 4; return Axis::X;
                case GateKind::RX:   angle = g.params[0]; return Axis::X;
                case GateKind::Y:    angle = kPi;      phase = kPi / 2;  return Axis::Y;  // Y = i RY(pi)
                case GateKind::RY:   angle = g.params[0]; return Axis::Y;
                default:             return Axis::None;
            }
        }

        struct NamedAngle { double angle; GateKind kind; double phase; };   // R(angle) = e^{i phase} kind

        constexpr NamedAngle kNamedZ[] = {
            {kPi, GateKind::Z, 0}, {kPi / 2, GateKind::S, 0}, {-kPi / 2, GateKind::Sdg, 0},
            {kPi / 4, GateKind::T, 0}, {-kPi / 4, GateKind::Tdg, 0},
        };
        constexpr NamedAngle kNamedX[] = {
            {kPi, GateKind::X, -kPi / 2}, {kPi / 2, GateKind::SX, -kPi / 4}, {-kPi / 2, GateKind::SXdg, kPi / 4},
        };
        constexpr NamedAngle kNamedY[] = {
            {kPi, GateKind::Y, -kPi / 2},
        };

        // Canonical gate for R(angle) about axis on qubit q, adding to phase so that
        // R(angle) = e^{i phase_added} out. Returns false if R(angle) is the identity.
        bool from_rotation(Axis axis, double angle, uint8_t q, Gate& out, double& phase) {
            const double k = wrap(angle);
            if (axis != Axis::Z) phase += k * kPi;    // RX/RY(t + 2pi) = -RX/RY(t); P has period 2pi
            if (near(angle, 0)) return false;

            auto pick = [&](const auto& table, GateKind fallback) {
                for (const NamedAngle& na : table)
                    if (near(angle, na.angle)) { out = make_gate(na.kind, q); phase += na.phase; return; }
                out = make_gate(fallback, q, {angle, 0, 0});
            };
            switch (axis) {
                case Axis::Z: pick(kNamedZ, GateKind::P);  break;
                case Axis::X: pick(kNamedX, GateKind::RX); break;
                case Axis::Y: pick(kNamedY, GateKind::RY); break;
                default:      return false;
            }
            return true;
        }

        // Rewrites node id into canonical form (steps 1-3), removing it if it is the identity.
        void canonicalize(CircuitDAG& dag, uint32_t id) {
            const Gate g = dag.nodes[id].gate;

            double angle, phase;
            const Axis axis = as_rotation(g, angle, phase);
            if (axis != Axis::None) {
                Gate out;
                if (from_rotation(axis, angle, g.q[0], out, phase)) dag.replace(id, out);
                else dag.remove(id);
                dag.global_phase += phase;
                return;
            }

            if (g.kind == GateKind::U) {
                double theta = g.params[0], phi = g.params[1], lam = g.params[2];
                dag.global_phase += wrap(theta) * kPi;                       // U(t + 2pi, ..) = -U(t, ..)
                if (theta < 0) { theta = -theta; phi += kPi; lam += kPi; }   // U(-t, f, l) = U(t, f + pi, l + pi)
                wrap(phi);
                wrap(lam);
                dag.replace(id, make_gate(GateKind::U, g.q[0], {theta, phi, lam}));
                return;
            }

            if (info(g.kind).symmetric && g.q[0] > g.q[1])
                dag.replace(id, make_gate(g.kind, g.q[1], g.q[0]));
        }

        // The node directly after id on every one of its wires, if it acts on exactly the same
        // qubits; kNone otherwise.
        uint32_t sole_successor(const CircuitDAG& dag, uint32_t id) {
            const CircuitDAG::Node& nd = dag.nodes[id];
            const uint32_t b = nd.next[0];
            if (b == kNone || dag.nodes[b].gate.nq != nd.gate.nq) return kNone;
            if (nd.gate.nq == 2 && nd.next[1] != b) return kNone;
            return b;
        }

        void push_preds(const CircuitDAG& dag, uint32_t id, std::vector<uint32_t>& work) {
            const CircuitDAG::Node& nd = dag.nodes[id];
            for (int i = 0; i < nd.gate.nq; ++i)
                if (nd.prev[i] != kNone) work.push_back(nd.prev[i]);
        }

        // Fuses or cancels node a with its sole successor (steps 4-5). Any node whose successor
        // may have changed is pushed onto work.
        void try_merge(CircuitDAG& dag, uint32_t a, std::vector<uint32_t>& work) {
            const uint32_t b = sole_successor(dag, a);
            if (b == kNone) return;
            const Gate ga = dag.nodes[a].gate, gb = dag.nodes[b].gate;

            double ta, pa, tb, pb;
            const Axis axis = as_rotation(ga, ta, pa);
            if (axis != Axis::None && axis == as_rotation(gb, tb, pb)) {
                double phase = pa + pb;
                Gate out;
                push_preds(dag, a, work);
                dag.remove(b);
                if (from_rotation(axis, ta + tb, ga.q[0], out, phase)) {
                    dag.replace(a, out);
                    work.push_back(a);
                } else {
                    dag.remove(a);
                }
                dag.global_phase += phase;
                return;
            }

            if (same_gate(gb, inverse(ga), kAngleTol)) {
                push_preds(dag, a, work);
                dag.remove(a);
                dag.remove(b);
            }
        }
    }

    void normalize(CircuitDAG& dag) {
        const uint32_t count = static_cast<uint32_t>(dag.nodes.size());
        for (uint32_t id = 0; id < count; ++id)
            if (dag.nodes[id].alive) canonicalize(dag, id);

        // Every successful merge removes a node, so the worklist drains in O(nodes) merges.
        std::vector<uint32_t> work;
        for (uint32_t id = count; id-- > 0;)
            if (dag.nodes[id].alive) work.push_back(id);
        while (!work.empty()) {
            const uint32_t a = work.back();
            work.pop_back();
            if (dag.nodes[a].alive) try_merge(dag, a, work);
        }

        wrap(dag.global_phase);
    }

}
