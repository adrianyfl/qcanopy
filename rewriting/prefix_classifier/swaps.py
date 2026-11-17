"""Elide explicit, unconditional SWAP gates from an immutable wire DAG."""
from dataclasses import dataclass
from typing import Tuple

from .model import CircuitDAG, InputError, digest, integer, param_json


@dataclass(frozen=True)
class SwapRemovalResult:
    """Core circuit plus its indispensable original-output-to-core mapping.

    Input wire ordering is unchanged. For little-endian state vectors,
    original[i] = core[output_index_to_core(i)], including global phase.
    """

    core_dag: CircuitDAG
    output_wire_to_core: Tuple[int, ...]
    removed_gate_ids: Tuple[str, ...]
    input_dag_fingerprint: str

    def output_index_to_core(self, index: int) -> int:
        """Translate one original output basis index; q0 is the least significant bit."""
        integer(index, "output index")
        if index >= (1 << self.core_dag.num_qubits):
            raise InputError("output index is outside the state vector")
        return sum(((index >> wire) & 1) << core
                   for wire, core in enumerate(self.output_wire_to_core))

    def to_dict(self):
        """Serialize a wrapper, not a bare qdag/1; retain mapping with the core."""
        return {
            "schema": "qswap-removal/1",
            "input_dag_fingerprint": self.input_dag_fingerprint,
            "core_dag": self.core_dag.to_dict(),
            "output_wire_to_core": list(self.output_wire_to_core),
            "removed_gate_ids": list(self.removed_gate_ids),
        }


def remove_swaps(dag: CircuitDAG) -> SwapRemovalResult:
    """Return a fresh DAG and final output permutation without modifying input.

    The reserved, case-sensitive op 'swap' means the ordinary two-qubit SWAP
    with no parameters. Other gates retain their exact name, parameters and
    ordered operands. No CX reconstruction, QASM conversion or matching occurs.
    The caller supplies a validated, unitary qdag/1 (no classical control).
    """
    mapping = list(range(dag.num_qubits))
    last = [None] * dag.num_qubits
    nodes, edges, removed = [], [], []
    # Traverse the original DAG only; rewritten wire dependencies may differ.
    for gid in dag.order:
        gate = dag.by_id[gid]
        if gate.op == "swap":
            if len(gate.qubits) != 2 or gate.params:
                raise InputError(f"gate {gid}: swap requires two qubits and no parameters")
            a, b = gate.qubits
            mapping[a], mapping[b] = mapping[b], mapping[a]
            removed.append(gid)
            continue
        qubits = [mapping[q] for q in gate.qubits]
        nodes.append({"id": gid, "op": gate.op, "qubits": qubits,
                      "params": [param_json(p) for p in gate.params]})
        for q in qubits:
            if last[q] is not None:
                edges.append({"source": last[q], "target": gid, "qubit": q})
            last[q] = gid
    core = CircuitDAG.from_dict({
        "schema": "qdag/1", "num_qubits": dag.num_qubits,
        "global_phase": dag.global_phase, "nodes": nodes, "edges": edges,
    })
    return SwapRemovalResult(core, tuple(mapping), tuple(removed), digest(dag.to_dict()))
