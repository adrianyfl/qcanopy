"""Import QCanopy circuit JSON v1 without angle snapping or normalization."""
import json
import math
from pathlib import Path

from .model import CircuitDAG, InputError, fields, integer


# Generic matrix gates require a matrix-aware matcher and are deliberately rejected.
GATES = {**{op: (1, 0) for op in
             ('id', 'h', 'x', 'y', 'z', 's', 'sdg', 't', 'tdg', 'sx', 'sxdg')},
         **{op: (1, 1) for op in ('rx', 'ry', 'rz', 'p')}, 'u': (1, 3),
         **{op: (2, 0) for op in ('cx', 'cz', 'swap')}}


def sequence_to_dag(n, gates, phase=0):
    """Build immediate per-wire edges; exact strings are allowed in this helper."""
    integer(n, 'n', 1)
    nodes, edges, last = [], [], [None] * n
    for i, gate in enumerate(gates):
        gid = 'g%d' % i
        qs = gate['qubits']
        if not isinstance(qs, list) or not qs or len(set(qs)) != len(qs):
            raise InputError('invalid ordered gate operands')
        for q in qs:
            integer(q, 'qubit')
            if q >= n:
                raise InputError('qubit out of range')
        nodes.append(dict(id=gid, op=gate['name'], qubits=qs, params=gate['params']))
        for q in qs:
            if last[q] is not None:
                edges.append(dict(source=last[q], target=gid, qubit=q))
            last[q] = gid
    return CircuitDAG.from_dict(dict(schema='qdag/1', num_qubits=n,
                                     global_phase=phase, nodes=nodes, edges=edges))


def parse_qcanopy_ir(data, *, allow_opaque_matrices=False):
    """Parse the named-gate subset; preserve numeric parameters exactly as supplied."""
    fields(data, {'schema_version', 'n', 'global_phase', 'gates'}, where='QCanopy IR')
    if type(data['schema_version']) is not int or data['schema_version'] != 1:
        raise InputError('expected QCanopy schema_version=1')
    integer(data['n'], 'n', 1)
    if data['n'] > 32:
        raise InputError('QCanopy v1 supports at most 32 qubits')
    if type(data['global_phase']) not in (int, float) or not math.isfinite(data['global_phase']):
        raise InputError('global_phase must be finite')
    if not isinstance(data['gates'], list):
        raise InputError('gates must be an array')
    gates=[]
    for g in data['gates']:
        if allow_opaque_matrices and g.get('name') in ('unitary_1q','unitary_2q'):
            # Opaque tokens preserve EVERY matrix entry, not a rounded hash or
            # inferred named gate. They cannot match named library rules.
            # Input unitary validity remains the qdag contract; native generator
            # validates its own matrix gates. No numerical equivalence oracle.
            fields(g, {'name','qubits','params','matrix'},where='opaque gate')
            nq=1 if g['name']=='unitary_1q' else 2
            if not isinstance(g['qubits'],list) or len(g['qubits'])!=nq or g['params']!=[]:
                raise InputError('invalid opaque operands/parameters')
            matrix=g['matrix']
            if not isinstance(matrix,list) or len(matrix)!=(1<<(2*nq)):
                raise InputError('invalid opaque matrix dimensions')
            params=[]
            for z in matrix:
                if not isinstance(z,list) or len(z)!=2 or any(type(v) not in (int,float) or not math.isfinite(v) for v in z):
                    raise InputError('opaque matrix entries must be finite [real,imaginary]')
                params.extend(z)
            gates.append(dict(name=g['name'],qubits=g['qubits'],params=params))
            continue
        fields(g, {'name', 'qubits', 'params'}, where='QCanopy gate')
        if not isinstance(g['name'], str) or g['name'] not in GATES:
            raise InputError('unsupported QCanopy gate: %s' % g['name'])
        nq, np = GATES[g['name']]
        if not isinstance(g['qubits'], list) or len(g['qubits']) != nq:
            raise InputError('wrong gate arity')
        if not isinstance(g['params'], list) or len(g['params']) != np:
            raise InputError('wrong parameter count')
        if any(type(p) not in (int, float) or not math.isfinite(p) for p in g['params']):
            raise InputError('QCanopy parameters must be finite numbers')
        gates.append(g)
    return sequence_to_dag(data['n'], gates, data['global_phase'])


def load_qcanopy_ir(path):
    return parse_qcanopy_ir(json.loads(Path(path).read_text()))
