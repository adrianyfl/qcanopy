"""Causally closed boundary windows; no artificial edges across removed gates."""
from .model import CircuitDAG, integer, param_json


def boundary_windows(dag, max_depth=24):
    """Return prefix/downset and suffix/upset DAGs, preserving IDs and wire labels.

    Depth is longest-path gate depth in the ORIGINAL full DAG. Every true prefix
    (suffix) region of <= max_depth gates is contained in the corresponding
    window: its full predecessor (successor) closure is the region itself.
    This is a boundary search scope, not an exhaustive search for middle matches.
    """
    integer(max_depth, 'max_depth', 1)
    head, tail = {}, {}
    for g in dag.order:
        head[g] = 1 + max((head[p] for p in dag.predecessors[g]), default=0)
    for g in reversed(dag.order):
        tail[g] = 1 + max((tail[p] for p in dag.successors[g]), default=0)
    result = {}
    for side, depth in (('prefix', head), ('suffix', tail)):
        ids = {g for g in dag.order if depth[g] <= max_depth}
        result[side] = CircuitDAG.from_dict(dict(schema='qdag/1', num_qubits=dag.num_qubits,
            global_phase=dag.global_phase,
            nodes=[dict(id=g, op=dag.by_id[g].op, qubits=list(dag.by_id[g].qubits),
                        params=[param_json(p) for p in dag.by_id[g].params]) for g in dag.order if g in ids],
            edges=[dict(source=a, target=b, qubit=q) for a,b,q in dag.edges if a in ids and b in ids]))
    return result
