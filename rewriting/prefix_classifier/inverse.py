"""Lazy identity schemas G G† = I, with exact parameter relations and DAG cuts.

Supported gates are declared in inverse_schemas.json. Gate-group cancellation
follows by repeated adjacent-pair elimination, not an expanding class catalog.
"""
import ast
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path
import re

from .model import CircuitDAG, Gate, canonical_json, digest, param_json


@lru_cache(maxsize=1)
def schemas():
    return json.loads(Path(__file__).with_name('inverse_schemas.json').read_text())


def angle_key(value):
    """Exact affine expression over real symbols, including pi. No eval/trig.

    Floats are their exact binary rationals; string decimals are exact decimals.
    Nonlinear expressions/calls are deliberately unsupported, not approximated.
    """
    if type(value) in (int, float):
        try:
            q=Fraction(value)
            return (('',q),) if q else ()
        except (ValueError,OverflowError):return None
    if not isinstance(value,str):return None
    source=value.strip()
    if len(source)>2048:return None
    try:
        tree=ast.parse(source,mode='eval')
        if sum(1 for _ in ast.walk(tree))>128:return None
        def scale(a,q):return {k:v*q for k,v in a.items() if v*q}
        def visit(n):
            if isinstance(n,ast.Constant) and type(n.value) in (int,float):
                token=ast.get_source_segment(source,n).replace('_','')
                exponent=re.search(r'[eE]([+-]?\d+)$',token)
                if len(token)>256 or (exponent and abs(int(exponent[1]))>4096):
                    raise ValueError('angle literal too large')
                q=Fraction(token)
                return {'':q} if q else {}
            if isinstance(n,ast.Name) and n.id not in ('nan','inf','Infinity'):
                return {n.id:Fraction(1)}
            if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd)):
                return scale(visit(n.operand),-1 if isinstance(n.op,ast.USub) else 1)
            if isinstance(n,ast.BinOp):
                a,b=visit(n.left),visit(n.right)
                if isinstance(n.op,(ast.Add,ast.Sub)):
                    out=dict(a);sign=-1 if isinstance(n.op,ast.Sub) else 1
                    for k,v in b.items():out[k]=out.get(k,0)+sign*v
                    return {k:v for k,v in out.items() if v}
                if isinstance(n.op,ast.Mult):
                    if set(a)<={''}:return scale(b,a.get('',0))
                    if set(b)<={''}:return scale(a,b.get('',0))
                if isinstance(n.op,ast.Div) and set(b)<={''} and b.get('',0):
                    return scale(a,1/b[''])
            raise ValueError('unsupported angle expression')
        return tuple(sorted(visit(tree.body).items()))
    except (ValueError,SyntaxError,TypeError,ZeroDivisionError,RecursionError):return None


def _negative(value):
    if angle_key(value) is None:return None
    return -value if type(value) in (int,float) else '-('+value+')'


def dagger_gate(gate):
    """Return a recognized exact dagger descriptor, or None for unsupported ops."""
    rules=schemas();fixed=rules['fixed'];angles=rules['negate_angle']
    if gate.op in fixed:
        op,n=fixed[gate.op]
        if len(gate.qubits)!=n or gate.params:return None
        return Gate(gate.id,op,gate.qubits,())
    if gate.op in angles:
        if len(gate.qubits)!=angles[gate.op] or len(gate.params)!=1:return None
        p=_negative(gate.params[0])
        return None if p is None else Gate(gate.id,gate.op,gate.qubits,(p,))
    if gate.op=='u' and len(gate.qubits)==1 and len(gate.params)==3:
        p=tuple(_negative(gate.params[i]) for i in (0,2,1))
        return None if any(x is None for x in p) else Gate(gate.id,'u',gate.qubits,p)
    return None


def inverse_relation(a,b):
    """An exact schema witness, retaining control/target order and full arity."""
    expected=dagger_gate(a)
    if expected is None or expected.op!=b.op or expected.qubits!=b.qubits or len(expected.params)!=len(b.params):
        return None
    if any(angle_key(p) is None or angle_key(p)!=angle_key(q) for p,q in zip(expected.params,b.params)):
        return None
    return dict(identity_class_id=schemas()['identity_class_id'],
                rule_id='dagger/'+a.op,lhs_gate_ids=[a.id,b.id],
                ordered_qubits=list(a.qubits),parameters_A=[param_json(p) for p in a.params],
                parameters_B=[param_json(p) for p in b.params])


@dataclass(frozen=True)
class InverseCancellationResult:
    core_dag: CircuitDAG
    witnesses: tuple
    input_fingerprint: str

    def to_dict(self, *, include_core=True):
        out=dict(schema='qinverse-cancellation/1',identity_class_id=schemas()['identity_class_id'],
            equivalence='exact_unitary',input_fingerprint=self.input_fingerprint,
            witnesses=list(self.witnesses),removed_gates=2*len(self.witnesses))
        if include_core:out['core_dag']=self.core_dag.to_dict()
        return out


def cancel_inverses(dag):
    """Cancel on all operand wires; disjoint gates may interleave safely.

    Single topological pass, O(gates + edges) stack operations. A witness list
    is a replayable sequence: each pair is adjacent after earlier pairs vanish.
    G1 G2 ... Gk Gk† ... G2† G1† reduces from the inside out. Unknown gates block
    their wires; no arbitrary matrix inverse or equivalence search is attempted.
    """
    active={};last=[None]*dag.num_qubits;witnesses=[]
    for gid in dag.order:
        gate=dag.by_id[gid];previous=last[gate.qubits[0]]
        witness=None
        if previous is not None and all(last[q]==previous for q in gate.qubits):
            witness=inverse_relation(active[previous][0],gate)
        if witness is not None:
            _,pred=active.pop(previous)
            for q in gate.qubits:last[q]=pred[q]
            witnesses.append(witness)
        else:
            active[gid]=(gate,{q:last[q] for q in gate.qubits})
            for q in gate.qubits:last[q]=gid
    if not witnesses:return InverseCancellationResult(dag,(),digest(dag.to_dict()))
    nodes=[];edges=[];last=[None]*dag.num_qubits
    for gid in dag.order:
        if gid not in active:continue
        g=dag.by_id[gid]
        nodes.append(dict(id=gid,op=g.op,qubits=list(g.qubits),params=[param_json(p) for p in g.params]))
        for q in g.qubits:
            if last[q] is not None:edges.append(dict(source=last[q],target=gid,qubit=q))
            last[q]=gid
    core=CircuitDAG.from_dict(dict(schema='qdag/1',num_qubits=dag.num_qubits,global_phase=dag.global_phase,nodes=nodes,edges=edges))
    return InverseCancellationResult(core,tuple(witnesses),digest(dag.to_dict()))


def replay_inverse_certificate(dag,certificate):
    """Independently check chronological adjacency and every exact schema step.

    Uses mutable wire links rather than the construction pass's stacks.
    Returns False for modified inputs, duplicate consumption or invalid cuts.
    """
    try:
        if certificate['identity_class_id']!=schemas()['identity_class_id'] or certificate['input_fingerprint']!=digest(dag.to_dict()):
            return False
        active=set(dag.order);prev=dict(dag.prev);nxt=dict(dag.next)
        for witness in certificate['witnesses']:
            a,b=witness['lhs_gate_ids']
            if a==b or a not in active or b not in active:return False
            ga,gb=dag.by_id[a],dag.by_id[b]
            if inverse_relation(ga,gb)!=witness:return False
            if any(nxt.get((a,q))!=b or prev.get((b,q))!=a for q in ga.qubits):return False
            for q in ga.qubits:
                left,right=prev.get((a,q)),nxt.get((b,q))
                if left is not None:
                    if right is None:nxt.pop((left,q),None)
                    else:nxt[left,q]=right
                if right is not None:
                    if left is None:prev.pop((right,q),None)
                    else:prev[right,q]=left
            active.remove(a);active.remove(b)
        if certificate['removed_gates']!=len(dag.nodes)-len(active):return False
        if 'core_dag' in certificate:
            core=CircuitDAG.from_dict(certificate['core_dag'])
            if core.num_qubits!=dag.num_qubits or set(core.by_id)!=active or core.global_phase!=dag.global_phase:return False
            for gid in active:
                if core.by_id[gid]!=dag.by_id[gid]:return False
            edges={(a,b,q) for (a,q),b in nxt.items() if a in active and b in active}
            if edges!=set(core.edges):return False
        return True
    except (KeyError,ValueError,TypeError):return False
