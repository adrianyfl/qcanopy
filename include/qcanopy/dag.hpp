#pragma once

#include <array>
#include <complex>
#include <cstdint>
#include <string_view>
#include <vector>
#include "qcanopy/ir.hpp"

namespace qcanopy {

    struct CircuitDAG {
        static constexpr uint32_t kNone = UINT32_MAX;

        struct Node {
            Gate gate;
            std::array<uint32_t, 2> prev{kNone, kNone};  // prev[i]: previous node on wire gate.q[i]
            std::array<uint32_t, 2> next{kNone, kNone};  // next[i]: next node on wire gate.q[i]
            bool alive = true;                           // tombstone; compact later
        };

        int n = 0;
        std::vector<Node> nodes;                         // one contiguous pool
        std::array<uint32_t, kMaxQubits> head{kNone};    // first node on each wire
        std::array<uint32_t, kMaxQubits> tail{kNone};    // last node on each wire
        double global_phase = 0.0;

        explicit CircuitDAG(const Circuit& c);
        Circuit to_circuit() const;                      // topological walk over alive nodes

        uint32_t next_on(uint32_t id, uint8_t q) const;  // O(1)
        uint32_t prev_on(uint32_t id, uint8_t q) const;  // O(1)
        void remove(uint32_t id);                        // splice out of both wires, O(1)
        void replace(uint32_t id, const Gate& g);        // same qubits, in place

        // Adding gates. Each returns the new node's id. For insert_after/insert_before,
        // every qubit of g must be one of anchor's qubits (this keeps the graph acyclic).
        uint32_t append(const Gate& g);                             // at the end of the circuit
        uint32_t insert_after(uint32_t anchor, const Gate& g);      // directly after anchor
        uint32_t insert_before(uint32_t anchor, const Gate& g);     // directly before anchor

    private:
        // Splice slot i of node id into wire gate.q[i], directly after node a (kNone = wire start).
        void link_after(uint32_t id, int i, uint32_t a);
    };

}

