"""Validated, immutable inputs for boundary classification (Python 3.9+)."""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import heapq
import json
import math
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping, Tuple


class InputError(ValueError):
    """An input violates the versioned contract."""


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      allow_nan=False)


def digest(value: Any) -> str:
    return sha256(canonical_json(value).encode()).hexdigest()


def fields(value, required, optional=(), where="object"):
    if not isinstance(value, dict):
        raise InputError(f"{where}: expected an object")
    missing = set(required) - value.keys()
    unknown = value.keys() - set(required) - set(optional)
    if missing or unknown:
        raise InputError(f"{where}: missing={sorted(missing)}, unknown={sorted(unknown)}")


def integer(value, where, minimum=0):
    if type(value) is not int or value < minimum:
        raise InputError(f"{where}: expected integer >= {minimum}")
    return value


def identifier(value, where):
    if not isinstance(value, str) or not value:
        raise InputError(f"{where}: expected a nonempty string")
    return value


def array(value, where):
    if not isinstance(value, list):
        raise InputError(f"{where}: expected an array")
    return value


@dataclass(frozen=True)
class Variable:
    name: str


def parameter(value, where, allow_variable=False):
    if isinstance(value, dict) and allow_variable:
        fields(value, {"var"}, where=where)
        return Variable(identifier(value["var"], where))
    if isinstance(value, str):
        return identifier(value, where)
    if type(value) in (int, float):
        if type(value) is float:
            if not math.isfinite(value):
                raise InputError(f"{where}: parameters must be finite")
            if value.is_integer():
                return int(value)
        return value
    raise InputError(f"{where}: expected a finite number, exact string, or declared variable")


def param_json(value):
    return {"var": value.name} if isinstance(value, Variable) else value


def param_equal(left, right):
    # No tolerance, symbolic evaluation, or angle rounding is performed.
    return type(left) is type(right) and left == right


@dataclass(frozen=True)
class Gate:
    id: str
    op: str
    qubits: Tuple[int, ...]
    params: Tuple[Any, ...] = ()


@dataclass(frozen=True)
class CircuitDAG:
    num_qubits: int
    nodes: Tuple[Gate, ...]
    edges: Tuple[Tuple[str, str, int], ...]
    global_phase: Any
    by_id: Mapping[str, Gate]
    prev: Mapping[Tuple[str, int], str]
    next: Mapping[Tuple[str, int], str]
    predecessors: Mapping[str, Tuple[str, ...]]
    successors: Mapping[str, Tuple[str, ...]]
    order: Tuple[str, ...]
    rank: Mapping[str, int]

    @classmethod
    def from_dict(cls, data, *, pattern=False):
        fields(data, {"schema", "num_qubits", "nodes", "edges"}, {"global_phase"}, "DAG")
        if data["schema"] != "qdag/1":
            raise InputError("DAG: expected schema qdag/1")
        n = integer(data["num_qubits"], "num_qubits", 1)
        nodes = []
        by_id = {}
        forbidden = {"measure", "reset", "barrier", "delay", "if_else", "while_loop",
                     "for_loop", "switch_case", "store", "kraus", "superop", "initialize"}
        for raw in array(data["nodes"], "nodes"):
            fields(raw, {"id", "op", "qubits"}, {"params"}, "gate")
            gid = identifier(raw["id"], "gate.id")
            op = identifier(raw["op"], "gate.op")
            if op.lower() in forbidden:
                raise InputError(f"gate {gid}: non-unitary/directive operation {op!r} is unsupported")
            qs = tuple(integer(q, f"gate {gid}.qubits") for q in array(raw["qubits"], "qubits"))
            if not qs or len(set(qs)) != len(qs) or any(q >= n for q in qs):
                raise InputError(f"gate {gid}: invalid ordered qubit operands")
            ps = tuple(parameter(p, f"gate {gid}.params", pattern)
                       for p in array(raw.get("params", []), "params"))
            if gid in by_id:
                raise InputError(f"duplicate gate id {gid!r}")
            gate = Gate(gid, op, qs, ps)
            nodes.append(gate)
            by_id[gid] = gate
        prev, nxt, edges = {}, {}, []
        pred = {g.id: set() for g in nodes}
        succ = {g.id: set() for g in nodes}
        for raw in array(data["edges"], "edges"):
            fields(raw, {"source", "target", "qubit"}, where="edge")
            a, b = raw["source"], raw["target"]
            identifier(a, "edge.source")
            identifier(b, "edge.target")
            q = integer(raw["qubit"], "edge.qubit")
            if a not in by_id or b not in by_id or a == b:
                raise InputError("edge: unknown endpoint or self edge")
            if q not in by_id[a].qubits or q not in by_id[b].qubits:
                raise InputError("edge: qubit must occur at both endpoints")
            if (a, q) in nxt or (b, q) in prev:
                raise InputError("edge: duplicate, branching, or merging wire edge")
            nxt[a, q], prev[b, q] = b, a
            succ[a].add(b)
            pred[b].add(a)
            edges.append((a, b, q))
        # Index once: rescanning all gates for every wire costs O(n_qubits * gates).
        used_by_wire = [set() for _ in range(n)]
        for gate in nodes:
            for q in gate.qubits:
                used_by_wire[q].add(gate.id)
        # A circuit wire is a single chain, not an arbitrary dependency graph.
        for q in range(n):
            used = used_by_wire[q]
            if not used:
                if pattern:
                    raise InputError("patterns must use every declared local qubit")
                continue
            heads = [g for g in used if (g, q) not in prev]
            if len(heads) != 1:
                raise InputError(f"qubit {q}: expected one complete wire chain")
            seen, node = set(), heads[0]
            while node is not None and node not in seen:
                seen.add(node)
                node = nxt.get((node, q))
            if node is not None or seen != used:
                raise InputError(f"qubit {q}: cycle or disconnected wire chain")
        indegree = {g: len(ps) for g, ps in pred.items()}
        ready = [g for g, degree in indegree.items() if degree == 0]
        heapq.heapify(ready)
        order = []
        while ready:
            a = heapq.heappop(ready)
            order.append(a)
            for b in sorted(succ[a]):
                indegree[b] -= 1
                if indegree[b] == 0:
                    heapq.heappush(ready, b)
        if len(order) != len(nodes):
            raise InputError("DAG contains a cycle across wires")
        if pattern and not nodes:
            raise InputError("empty pattern members are unsupported in qclass-library/1")
        phase = parameter(data.get("global_phase", 0), "global_phase")
        if pattern and phase != 0:
            raise InputError("v1 patterns must have zero explicit global_phase")
        return cls(n, tuple(sorted(nodes, key=lambda g: g.id)), tuple(sorted(edges)), phase,
                   MappingProxyType(by_id), MappingProxyType(prev), MappingProxyType(nxt),
                   MappingProxyType({g: tuple(sorted(v)) for g, v in pred.items()}),
                   MappingProxyType({g: tuple(sorted(v)) for g, v in succ.items()}),
                   tuple(order), MappingProxyType({g: i for i, g in enumerate(order)}))

    def to_dict(self):
        return {"schema": "qdag/1", "num_qubits": self.num_qubits,
                "global_phase": self.global_phase,
                "nodes": [{"id": g.id, "op": g.op, "qubits": list(g.qubits),
                           "params": [param_json(p) for p in g.params]} for g in self.nodes],
                "edges": [{"source": a, "target": b, "qubit": q} for a, b, q in self.edges]}

    def boundary(self, gates):
        gates = frozenset(gates)
        return (all(p in gates for g in gates for p in self.predecessors[g]),
                all(s in gates for g in gates for s in self.successors[g]))

    def ports(self, gates):
        gates = frozenset(gates)
        first, last = {}, {}
        for gid in sorted(gates, key=self.rank.__getitem__):
            for q in self.by_id[gid].qubits:
                first.setdefault(q, gid)
                last[q] = gid
        return first, last


@dataclass(frozen=True)
class Member:
    id: str
    dag: CircuitDAG


@dataclass(frozen=True)
class EquivalenceClass:
    id: str
    num_qubits: int
    parameters: Tuple[str, ...]
    members: Tuple[Member, ...]


@dataclass(frozen=True)
class ClassLibrary:
    classes: Tuple[EquivalenceClass, ...]
    equivalence: str = "exact_unitary"

    @classmethod
    def from_dict(cls, data):
        fields(data, {"schema", "equivalence", "classes"}, where="class library")
        if data["schema"] != "qclass-library/1" or data["equivalence"] not in (
                "exact_unitary", "up_to_global_phase"):
            raise InputError("expected qclass-library/1 with exact_unitary or up_to_global_phase equivalence")
        classes, ids = [], set()
        for raw in array(data["classes"], "classes"):
            fields(raw, {"id", "num_qubits", "members"}, {"parameters"}, "class")
            cid = identifier(raw["id"], "class.id")
            if cid in ids:
                raise InputError(f"duplicate class id {cid!r}")
            ids.add(cid)
            n = integer(raw["num_qubits"], "class.num_qubits", 1)
            parameters = tuple(identifier(p, "parameter name")
                               for p in array(raw.get("parameters", []), "parameters"))
            if len(set(parameters)) != len(parameters):
                raise InputError("duplicate parameter declaration")
            members, mids = [], set()
            for member in array(raw["members"], "members"):
                fields(member, {"id", "dag"}, where="member")
                mid = identifier(member["id"], "member.id")
                if mid in mids:
                    raise InputError(f"duplicate member id {mid!r} in {cid}")
                mids.add(mid)
                dag = CircuitDAG.from_dict(member["dag"], pattern=True)
                if dag.num_qubits != n:
                    raise InputError(f"member {mid}: arity differs from class")
                used = {p.name for g in dag.nodes for p in g.params if isinstance(p, Variable)}
                if used != set(parameters):
                    raise InputError(f"member {mid}: variables must match class parameter declarations")
                members.append(Member(mid, dag))
            if not members:
                raise InputError(f"class {cid}: expected at least one supplied member")
            classes.append(EquivalenceClass(cid, n, parameters,
                                            tuple(sorted(members, key=lambda m: m.id))))
        return cls(tuple(sorted(classes, key=lambda c: c.id)), data["equivalence"])

    def to_dict(self):
        return {"schema": "qclass-library/1", "equivalence": self.equivalence, "classes": [
            {"id": c.id, "num_qubits": c.num_qubits, "parameters": list(c.parameters),
             "members": [{"id": m.id, "dag": m.dag.to_dict()} for m in c.members]}
            for c in self.classes]}


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise InputError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def read_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream, object_pairs_hook=_unique_object)


def load_dag(path):
    return CircuitDAG.from_dict(read_json(path))


def load_library(path):
    return ClassLibrary.from_dict(read_json(path))
