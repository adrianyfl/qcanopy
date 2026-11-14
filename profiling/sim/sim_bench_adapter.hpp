#pragma once
// The only place sim_bench touches QCanopy APIs beyond sim.hpp.
//
// Everything sim_bench needs from the rest of the repository is these four
// calls. If an upstream signature differs from what is assumed here, fix it
// in this file only.
//
//   read_ir           qcanopy::testgen::read_ir(std::istream&) -> Circuit
//                     (documented in docs/generator.md, "Library use")
//   normalized_copy   qcanopy::testgen::normalized_copy(const Circuit&) -> Circuit
//   make_dag          CircuitDAG constructed from a Circuit
//   gate fields       Gate::nq, Gate::q[0..1], Gate::m (row-major, 4 or 16 entries),
//                     Circuit::n, Circuit::gates, Circuit::global_phase
//                     (exactly what sim.cpp already reads)

#include <istream>
#include <utility>

#include "qcanopy/dag.hpp"
#include "qcanopy/dataset_io.hpp"
#include "qcanopy/generator.hpp"
#include "qcanopy/ir.hpp"

namespace qcanopy::simbench {

inline Circuit load_ir(std::istream& in) { return qcanopy::testgen::read_ir(in); }

inline Circuit normalized(const Circuit& c) { return qcanopy::testgen::normalized_copy(c); }

inline CircuitDAG make_dag(const Circuit& c) { return CircuitDAG(c); }

}  // namespace qcanopy::simbench
