#include "qcanopy/dag.hpp"

#include<assert.h>

namespace qcanopy {

    CircuitDAG::CircuitDAG(const Circuit& c) : n(c.n), global_phase(c.global_phase) {
        head.fill(kNone);
        tail.fill(kNone);
        nodes.reserve(c.size());
        for (const Gate& g : c.gates) append(g);
    }

    Circuit CircuitDAG::to_circuit() const {
        return {}; // TODO
    }

    // find which slot of g holds qubit q. q must be one of g's qubits.
    inline int slot_of(const Gate& g, uint8_t q) {
        assert(g.q[0] == q || (g.nq == 2 && g.q[1] == q));
        return g.q[0] == q ? 0 : 1;
    }

    // return id of the next gate
    uint32_t CircuitDAG::next_on(uint32_t id, uint8_t q) const {
        const Node& nd = nodes[id];
        return nd.next[slot_of(nd.gate, q)];
    }

    uint32_t CircuitDAG::prev_on(uint32_t id, uint8_t q) const {
        const Node& nd = nodes[id];
        return nd.prev[slot_of(nd.gate, q)];
    }

    void CircuitDAG::link_after(uint32_t id, int i, uint32_t a) {
        uint8_t q = nodes[id].gate.q[i];
        uint32_t b = (a == kNone) ? head[q] : next_on(a, q);   // node that will follow us

        nodes[id].prev[i] = a;
        nodes[id].next[i] = b;

        if (a == kNone) head[q] = id;
        else            nodes[a].next[slot_of(nodes[a].gate, q)] = id;

        if (b == kNone) tail[q] = id;
        else            nodes[b].prev[slot_of(nodes[b].gate, q)] = id;
    }

    // Note: push_back may reallocate `nodes`, so never hold a Node& across it; use ids.
    uint32_t CircuitDAG::append(const Gate& g) {
        uint32_t id = static_cast<uint32_t>(nodes.size());
        nodes.push_back(Node{g});
        for (int i = 0; i < g.nq; ++i)
            link_after(id, i, tail[g.q[i]]);
        return id;
    }

    uint32_t CircuitDAG::insert_after(uint32_t anchor, const Gate& g) {
        assert(nodes[anchor].alive);
        uint32_t id = static_cast<uint32_t>(nodes.size());
        nodes.push_back(Node{g});
        for (int i = 0; i < g.nq; ++i)
            link_after(id, i, anchor);
        return id;
    }

    uint32_t CircuitDAG::insert_before(uint32_t anchor, const Gate& g) {
        assert(nodes[anchor].alive);
        uint32_t id = static_cast<uint32_t>(nodes.size());
        nodes.push_back(Node{g});
        for (int i = 0; i < g.nq; ++i)
            link_after(id, i, prev_on(anchor, g.q[i]));
        return id;
    }

    void CircuitDAG::remove(uint32_t id) {
        Node& nd = nodes[id];
        assert(nd.alive);

        for (int i =0; i < nd.gate.nq; i ++) {
            uint8_t q = nd.gate.q[i];
            uint32_t a = nd.prev[i];
            uint32_t b = nd.next[i];

            if (a == kNone) head[q] = b;
            else nodes[a].next[slot_of(nodes[a].gate, q)] = b;

            if (b == kNone) tail[q] = a;
            else nodes[b].prev[slot_of(nodes[b].gate, q)] = a;
        }

        // mark as dead
        nd.alive = false;
    }

    void CircuitDAG::replace(uint32_t id, const Gate& g) {
        Node& nd = nodes[id];
        const Gate& old = nd.gate;
        assert(nd.alive);
        assert(g.nq == old.nq);

        if (g.nq == 2 && g.q[0] == old.q[1] && g.q[1] == old.q[0]) {
            std::swap(nd.prev[0], nd.prev[1]);
            std::swap(nd.next[0], nd.next[1]);
        } else {
            assert(g.q[0] == old.q[0] && (g.nq == 1 || g.q[1] == old.q[1]));
        }

        nd.gate = g;
    }

} 
