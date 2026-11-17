"""Sound, bounded pair proof by boundary peeling and the existing fixed library.

Only returns equivalent_up_to_global_phase or inconclusive. Never uses generator
labels/proofs, simulation, dense matrices, or numeric distance as a proof oracle.
"""
from collections import deque
from dataclasses import asdict, replace
from types import SimpleNamespace
import time

from .model import CircuitDAG, canonical_json, param_equal
from .engine import Limits, classify
from .overlap import OverlapLimits, search_overlaps, _replace
from .inverse import cancel_inverses


class _Frontier:
    def __init__(self, dag):
        self.dag=dag;self.active=set(dag.order)
        self.head=[None]*dag.num_qubits;self.tail=[None]*dag.num_qubits
        for gid in dag.order:
            for q in dag.by_id[gid].qubits:
                if self.head[q] is None:self.head[q]=gid
                self.tail[q]=gid

    def remove(self, ids):
        self.active.difference_update(ids)
        for q in range(self.dag.num_qubits):
            while self.head[q] is not None and self.head[q] not in self.active:
                self.head[q]=self.dag.next.get((self.head[q],q))
            while self.tail[q] is not None and self.tail[q] not in self.active:
                self.tail[q]=self.dag.prev.get((self.tail[q],q))

    def snapshot(self):
        # Only boundary regions have been removed. Every remaining wire is a
        # contiguous interval, so filtering original edges creates no shortcuts.
        data=self.dag.to_dict()
        data['nodes']=[g for g in data['nodes'] if g['id'] in self.active]
        data['edges']=[e for e in data['edges'] if e['source'] in self.active and e['target'] in self.active]
        return CircuitDAG.from_dict(data)


def _same(a,b):
    return a.op==b.op and a.qubits==b.qubits and len(a.params)==len(b.params) and all(param_equal(x,y) for x,y in zip(a.params,b.params))


def _peel(a,b,trace):
    removed=0
    for side in ('prefix','suffix'):
        ends_a=a.head if side=='prefix' else a.tail
        ends_b=b.head if side=='prefix' else b.tail
        todo=deque(range(a.dag.num_qubits));pairs=[]
        while todo:
            q=todo.popleft();x,y=ends_a[q],ends_b[q]
            if x is None or y is None:continue
            ga,gb=a.dag.by_id[x],b.dag.by_id[y]
            if any(ends_a[w]!=x for w in ga.qubits) or any(ends_b[w]!=y for w in gb.qubits):continue
            if not _same(ga,gb):continue
            a.remove([x]);b.remove([y]);pairs.append([x,y]);todo.extend(ga.qubits);removed+=2
        if pairs:trace.append(dict(kind='structural',side=side,gate_pairs=pairs))
    return removed


def structural_peel(left,right):
    if left.num_qubits!=right.num_qubits:raise ValueError('pair qubit counts differ')
    a,b=_Frontier(left),_Frontier(right);trace=[]
    _peel(a,b,trace)
    return dict(left=a.snapshot(),right=b.snapshot(),trace=trace,
                removed=len(left.nodes)+len(right.nodes)-len(a.active)-len(b.active),
                equivalent=not a.active and not b.active)


def _identity_classes(library):
    return {c.id for c in library.classes if any(all(g.op=='id' and not g.params and len(g.qubits)==1
        for g in m.dag.nodes) for m in c.members)}


def _certificates(result,side):
    sigs={g['key']:g['signature'] for g in result['composition_groups']}
    matches={m['id']:m for m in result['matches']}
    out=[]
    for r in result['regions']:
        if not r['boundary'][side]:continue
        signature=sigs[r['composition_key']]
        # The ordered physical ports are essential: CX(0,1) != CX(1,0).
        key=canonical_json(dict(signature=signature,physical_ports=r['qubit_map']))
        out.append((key,dict(region=r,signature=signature,matches=[matches[p] for p in r['parts']])))
    return out


def _semantic_step(a,b,library,limits,trace,statistics,deadline):
    dags=[a.snapshot(),b.snapshot()];results=[]
    for label,d in zip(('A','B'),dags):
        if not d.nodes:results.append(None);continue
        left=deadline-time.perf_counter()
        if left<=0:return False
        lim=replace(limits,time_limit_seconds=min(limits.time_limit_seconds,left))
        r=classify(d,library,lim);results.append(r)
        statistics.append(dict(side=label,stop=r['stop_reason'],statistics=r['statistics']))
    identities=_identity_classes(library)
    for side in ('prefix','suffix'):
        certs=[_certificates(r,side) if r else [] for r in results]
        for idx in (0,1):
            for _,cert in sorted(certs[idx],key=lambda x:-len(x[1]['region']['gate_ids'])):
                if all(atom['class_id'] in identities for atom in cert['signature']['atoms']):
                    (a,b)[idx].remove(cert['region']['gate_ids'])
                    trace.append(dict(kind='library_identity',side=side,circuit=('A','B')[idx],certificate=cert))
                    return True
        index={key:cert for key,cert in certs[1]}
        for key,ca in sorted(certs[0],key=lambda x:-len(x[1]['region']['gate_ids'])):
            cb=index.get(key)
            if cb is None:continue
            a.remove(ca['region']['gate_ids']);b.remove(cb['region']['gate_ids'])
            trace.append(dict(kind='library_pair',side=side,left=ca,right=cb))
            return True
    return False


def _inverse_pass(left,right):
    dags=[];proof=[]
    for label,d in zip(('A','B'),(left,right)):
        result=cancel_inverses(d)
        dags.append(result.core_dag)
        if result.witnesses:
            proof.append(dict(kind='inverse_identity',circuit=label,
                              certificate=result.to_dict(include_core=False)))
    return dags[0],dags[1],proof


def _base_proof(left,right,library,limits,deadline,max_rounds=32,enable_inverse=True):
    a,b=_Frontier(left),_Frontier(right);trace=[];statistics=[]
    for _ in range(max_rounds):
        _peel(a,b,trace)
        if not a.active and not b.active:return True,a.snapshot(),b.snapshot(),trace,statistics
        if time.perf_counter()>=deadline:break
        if enable_inverse:
            na,nb,extra=_inverse_pass(a.snapshot(),b.snapshot())
            if extra:
                a,b=_Frontier(na),_Frontier(nb);trace.extend(extra)
                continue
        if not _semantic_step(a,b,library,limits,trace,statistics,deadline):break
    return False,a.snapshot(),b.snapshot(),trace,statistics


def prove_pair(left,right,library,limits=None,overlap_limits=None,*,enable_overlap=False,
               base_seconds=6,overlap_seconds=8,max_alternatives=32,max_residual_gates=128,on_base=None,
               enable_inverse=True):
    """Proof-first comparison on fixed physical wires, allowing global phase.

    Pair extraction is conservative: equal class/composition payload AND equal
    ordered physical ports. Unresolved pairs are never labeled non-equivalent.
    ON explores existing proof-backed overlap alternatives, then retries OFF.
    """
    if left.num_qubits!=right.num_qubits:raise ValueError('pair qubit counts differ')
    if type(enable_inverse) is not bool:raise ValueError('enable_inverse must be boolean')
    limits=limits or Limits(max_qubits=3,max_gates=16,max_matches=10000,max_regions=20000,
                           max_match_attempts=1000000,time_limit_seconds=2)
    overlap_limits=overlap_limits or OverlapLimits(max_qubits=3,max_original_gates=12,
        max_virtual_gates=16,max_members_per_class=8,max_views_per_match=24,max_views=3000,
        max_candidates=1000,max_match_attempts=1000000,time_limit_seconds=2)
    started=time.perf_counter();plain=structural_peel(left,right)
    if enable_inverse:
        # Linear cancellation runs before the residual-size search cap, so a
        # long G G† need not be enumerated as a large library pattern.
        plain['left'],plain['right'],extra=_inverse_pass(plain['left'],plain['right'])
        plain['trace'].extend(extra)
    if len(plain['left'].nodes)+len(plain['right'].nodes)>max_residual_gates:
        row=dict(verdict='inconclusive',reason='residual_size_limit',structural_removed=plain['removed'],
                    residual_left=plain['left'].to_dict(),residual_right=plain['right'].to_dict(),proof=plain['trace'],
                    base_seconds=time.perf_counter()-started,overlap_seconds=0,base_proved=False,overlap_attempted=False)
        if on_base:on_base(dict(row))
        return row
    deadline=started+base_seconds
    ok,a,b,proof,stats=_base_proof(plain['left'],plain['right'],library,limits,deadline,enable_inverse=enable_inverse)
    proof=plain['trace']+proof
    base_time=time.perf_counter()-started
    row=dict(verdict='equivalent_up_to_global_phase' if ok else 'inconclusive',base_proved=ok,
        structural_removed=plain['removed'],base_residual_gates=[len(a.nodes),len(b.nodes)],
        base_seconds=base_time,overlap_seconds=0,base_statistics=stats,proof=proof,
        residual_left=a.to_dict(),residual_right=b.to_dict(),overlap_attempted=False)
    if on_base:on_base(dict(row))
    if ok or not enable_overlap:return row
    ov_start=time.perf_counter();deadline=ov_start+overlap_seconds;metadata=[];examined=0
    for index,d in enumerate((a,b)):
        if not d.nodes:continue
        left_seconds=deadline-time.perf_counter()
        if left_seconds<=0:break
        base=classify(d,library,replace(limits,time_limit_seconds=min(limits.time_limit_seconds,left_seconds)))
        remaining=deadline-time.perf_counter()
        if remaining<=0:break
        meta,candidates=search_overlaps(d,library,base,replace(overlap_limits,time_limit_seconds=min(overlap_limits.time_limit_seconds,remaining)))
        metadata.append(dict(circuit=('A','B')[index],search=meta))
        # Prefer fewer gates, then deterministic structure. This uses no source
        # transform, expected label, or intended representative from a manifest.
        candidates=sorted(candidates,key=lambda c:(len(c['alternative_dag']['nodes']),c['alternative_key']))
        for c in candidates:
            if examined>=max_alternatives or time.perf_counter()>=deadline:break
            if not any(c['boundary'].values()):continue
            examined+=1
            alt=CircuitDAG.from_dict(c['alternative_dag'])
            changed,_=_replace(d,c['original_gate_ids'],SimpleNamespace(dag=alt),tuple(c['qubit_map']),{},'@pair:')
            na,nb=(changed,b) if index==0 else (a,changed)
            success,ra,rb,newproof,extra=_base_proof(na,nb,library,limits,deadline,enable_inverse=enable_inverse)
            if success:
                row.update(verdict='equivalent_up_to_global_phase',proof=proof+[
                    dict(kind='overlap_substitution',circuit=('A','B')[index],candidate=c)]+newproof,
                    residual_left=ra.to_dict(),residual_right=rb.to_dict())
                ok=True;break
        if ok or examined>=max_alternatives:break
    row.update(overlap_attempted=True,overlap_seconds=time.perf_counter()-ov_start,
               overlap_statistics=metadata,alternatives_examined=examined)
    return row
