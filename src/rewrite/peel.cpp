#include "qcanopy/peel.hpp"

#include <cmath>
#include <vector>
#include <stdexcept>

namespace qcanopy {

    namespace {
        constexpr uint32_t kNone = CircuitDAG::kNone;
        static constexpr double tol = 1e-9;

        // match gates up to a phase, true if b = e^{i alpha} a on the same qubits
        bool match_up_to_phase(const Gate& a, const Gate& b, double& alpha) {
            // check if qubits don't match
            if (a.nq != b.nq || a.q[0] != b.q[0] || (a.nq == 2 && a.q[1] != b.q[1])) return false;
            
            const int d = a.nq == 1 ? 4 : 16;
            int k = 0;
            for (int i = 1; i < d; ++i)
                if (std::abs(a.m[i]) > std::abs(a.m[k])) k = i;
            const cplx r = b.m[k] / a.m[k];
            if (std::abs(std::abs(r) - 1) > tol) return false;
            for (int i = 0; i < d; ++i)
                if (std::abs(b.m[i] - r * a.m[i]) > tol) return false;
            alpha = std::arg(r);
            return true;
        }
        
        // front or end node of wire q
        uint32_t boundary_node(const CircuitDAG& dag, uint8_t q, bool front) {
            const auto& ends = front ? dag.head : dag.tail;
            const uint32_t id = ends[q];
            if (id == kNone) return kNone;
            const Gate& g = dag.nodes[id].gate;
            for (int i = 0; i < g.nq; ++i)
                if (ends[g.q[i]] != id) return kNone;
            return id;
        }

        // peel matching gates from the front or end, returns the number of pairs removed
        std::size_t peel_boundary(CircuitDAG& R, CircuitDAG& C, bool front) {
            std::size_t removed = 0;
            std::vector<uint8_t> work;
            for (int q = R.n; q-- > 0;) work.push_back(static_cast<uint8_t>(q));

            while (!work.empty()) {
                const uint8_t q = work.back();
                work.pop_back();
                const uint32_t r = boundary_node(R, q, front);
                const uint32_t c = boundary_node(C, q, front);
                if (r == kNone || c == kNone) continue;

                double alpha;
                const Gate g = R.nodes[r].gate;
                if (!match_up_to_phase(g, C.nodes[c].gate, alpha)) continue;
                R.remove(r);
                C.remove(c);
                C.global_phase += alpha;
                ++removed;
                for (int i = 0; i < g.nq; ++i) work.push_back(g.q[i]);
            }
            return removed;
        }
    }

    std::size_t peel(CircuitDAG& R, CircuitDAG& C) {
        if (R.n != C.n) throw std::invalid_argument("peel: circuits have different qubit counts");

        const std::size_t prefix = peel_boundary(R, C, true);
        const std::size_t suffix = peel_boundary(R, C, false);
        return prefix + suffix;
    }

}