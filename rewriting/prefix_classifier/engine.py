"""Fixed-library matching and occurrence-backed composition links.

No rewrite, class union, new class/member, or representative extraction occurs.
Only occurrences, link witnesses, and boundary-region records are accumulated.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass
import heapq
from itertools import permutations
import time

from .model import (CircuitDAG, ClassLibrary, InputError, Variable, canonical_json,
                    digest, param_equal)


@dataclass(frozen=True)
class Limits:
    max_qubits: int = 4
    max_gates: int = 64
    max_rounds: int = 64
    max_matches: int = 10000
    max_regions: int = 20000
    max_link_witnesses: int = 50000
    max_match_attempts: int = 1000000
    max_extension_attempts: int = 1000000
    time_limit_seconds: float = 5.0

    def validate(self):
        for key, value in asdict(self).items():
            if key == "time_limit_seconds":
                if type(value) not in (int, float) or not (0 < value < float("inf")):
                    raise InputError("time_limit_seconds must be positive and finite")
            elif type(value) is not int or value < (0 if key == "max_rounds" else 1):
                raise InputError(f"{key}: invalid limit")
        if self.max_qubits > 6:
            raise InputError("v1 caps max_qubits at 6 to bound wire canonicalization")
        if self.max_gates > 256:
            raise InputError("v1 caps max_gates at 256")


class _Stopped(Exception):
    pass


@dataclass(frozen=True)
class _Occurrence:
    id: str
    class_id: str
    member_id: str
    gates: frozenset
    qubits: tuple
    bindings: tuple
    prefix: bool
    suffix: bool


@dataclass(frozen=True)
class _Region:
    id: str
    parts: tuple
    gates: frozenset
    qubits: tuple
    signature: str
    prefix: bool
    suffix: bool
    round: int
    parent: object
    added: object


class _Run:
    def __init__(self, dag, library, limits, on_round):
        self.dag, self.library, self.limits = dag, library, limits
        self.on_round = on_round
        self.start = time.perf_counter()
        self.counts = {"match_attempts": 0, "extension_attempts": 0,
                       "link_witnesses": 0, "width_pruned": 0,
                       "size_pruned": 0, "members_outside_bounds": 0,
                       "classes_started": 0, "members_started": 0, "members_completed": 0}
        self.occurrences, self.regions, self.links = [], [], []
        self.occ_keys, self.region_keys, self.link_keys = {}, {}, {}
        self.signatures = {}
        self.rounds = []
        self.ports_cache = {}
        self.starts, self.ends = defaultdict(list), defaultdict(list)
        # Visit circuit positions from the two ends, independent of input array order.
        scan = []
        left, right = 0, len(dag.order) - 1
        while left <= right:
            scan.append(dag.order[left])
            if right != left:
                scan.append(dag.order[right])
            left, right = left + 1, right - 1
        self.gate_index = defaultdict(list)
        for gid in scan:
            gate = dag.by_id[gid]
            self.gate_index[gate.op, len(gate.qubits), len(gate.params)].append(gid)

    def check(self, counter=None):
        if time.perf_counter() - self.start >= self.limits.time_limit_seconds:
            raise _Stopped("time_limit")
        if counter:
            limit = getattr(self.limits, "max_" + counter)
            if self.counts[counter] >= limit:
                raise _Stopped(counter + "_limit")
            self.counts[counter] += 1

    def ports(self, gates):
        if gates not in self.ports_cache:
            self.ports_cache[gates] = self.dag.ports(gates)
        return self.ports_cache[gates]

    def convex(self, gates):
        # Wire contiguity alone misses paths that leave on one wire and return
        # through another. Reject every outside path between two member gates.
        ceiling = max(self.dag.rank[g] for g in gates)
        queue = [s for g in gates for s in self.dag.successors[g] if s not in gates]
        seen = set()
        while queue:
            self.check()
            node = queue.pop()
            if node in gates:
                return False
            if node in seen or self.dag.rank[node] > ceiling:
                continue
            seen.add(node)
            queue.extend(self.dag.successors[node])
        return True

    def match_member(self, eqclass, member):
        pattern = member.dag
        if pattern.num_qubits > self.limits.max_qubits or len(pattern.nodes) > self.limits.max_gates:
            self.counts["members_outside_bounds"] += 1
            return

        def visit(position, matched, wires, reverse, bindings):
            self.check()
            if position == len(pattern.order):
                gates = frozenset(matched.values())
                if self.convex(gates):
                    yield gates, tuple(wires[q] for q in range(pattern.num_qubits)), tuple(
                        (p, bindings[p]) for p in eqclass.parameters)
                return
            pid = pattern.order[position]
            p = pattern.by_id[pid]
            forced = set()
            for q in p.qubits:
                previous = pattern.prev.get((pid, q))
                if previous is not None:
                    successor = self.dag.next.get((matched[previous], wires[q]))
                    if successor is None:
                        return
                    forced.add(successor)
            if len(forced) > 1:
                return
            candidates = sorted(forced) if forced else self.gate_index[p.op, len(p.qubits), len(p.params)]
            for gid in candidates:
                self.check("match_attempts")
                if gid in matched.values():
                    continue
                gate = self.dag.by_id[gid]
                if (gate.op, len(gate.qubits), len(gate.params)) != (p.op, len(p.qubits), len(p.params)):
                    continue
                w, rev, params = dict(wires), dict(reverse), dict(bindings)
                valid = True
                for local, actual in zip(p.qubits, gate.qubits):
                    if (local in w and w[local] != actual) or (actual in rev and rev[actual] != local):
                        valid = False
                        break
                    w[local], rev[actual] = actual, local
                    previous = pattern.prev.get((pid, local))
                    if previous is not None and self.dag.prev.get((gid, actual)) != matched[previous]:
                        valid = False
                        break
                if not valid:
                    continue
                for expected, actual in zip(p.params, gate.params):
                    if isinstance(expected, Variable):
                        if expected.name in params and not param_equal(params[expected.name], actual):
                            valid = False
                            break
                        params[expected.name] = actual
                    elif not param_equal(expected, actual):
                        valid = False
                        break
                if valid:
                    yield from visit(position + 1, {**matched, pid: gid}, w, rev, params)

        yield from visit(0, {}, {}, {}, {})

    def seed(self):
        for eqclass in self.library.classes:
            self.counts["classes_started"] += 1
            for member in eqclass.members:
                self.check()
                self.counts["members_started"] += 1
                for gates, wires, params in self.match_member(eqclass, member):
                    key = (eqclass.id, gates, wires, params)
                    if key in self.occ_keys:
                        continue
                    if len(self.occurrences) >= self.limits.max_matches:
                        raise _Stopped("matches_limit")
                    prefix, suffix = self.dag.boundary(gates)
                    occurrence = _Occurrence(f"m{len(self.occurrences)}", eqclass.id, member.id,
                                             gates, wires, params, prefix, suffix)
                    self.occ_keys[key] = len(self.occurrences)
                    self.occurrences.append(occurrence)
                    first, last = self.ports(gates)
                    for q, g in first.items():
                        self.starts[g, q].append(len(self.occurrences) - 1)
                    for q, g in last.items():
                        self.ends[g, q].append(len(self.occurrences) - 1)
                    # Interior seeds are needed to bridge forks/joins. Their
                    # boundary flags remain false until all dependencies close.
                    self.register((len(self.occurrences) - 1,), 0)
                self.counts["members_completed"] += 1

    def signature(self, parts):
        atoms = [self.occurrences[i] for i in parts]
        qubits = sorted({q for atom in atoms for q in atom.qubits})
        predecessors = [set() for _ in atoms]
        successors = [set() for _ in atoms]
        last = {}
        for i, atom in enumerate(atoms):
            for q in atom.qubits:
                if q in last:
                    predecessors[i].add(last[q])
                    successors[last[q]].add(i)
                last[q] = i
        best, best_wires, best_payload = None, None, None
        # Bounded wire permutations canonicalize renaming and independently
        # ordered atoms, not alternatives inside equivalence classes.
        for order in permutations(qubits):
            self.check()
            mapping = {q: i for i, q in enumerate(order)}
            payloads = [{"class_id": a.class_id, "bindings": dict(a.bindings),
                         "qubits": [mapping[q] for q in a.qubits]} for a in atoms]
            labels = [canonical_json(p) for p in payloads]
            degrees = [len(ps) for ps in predecessors]
            ready = [(labels[i], i) for i, degree in enumerate(degrees) if degree == 0]
            heapq.heapify(ready)
            sequence = []
            while ready:
                _, i = heapq.heappop(ready)
                sequence.append(payloads[i])
                for j in successors[i]:
                    degrees[j] -= 1
                    if not degrees[j]:
                        heapq.heappush(ready, (labels[j], j))
            payload = {"num_qubits": len(order), "atoms": sequence}
            encoded = canonical_json(payload)
            if best is None or encoded < best:
                best, best_wires, best_payload = encoded, tuple(order), payload
        return digest(best_payload), best_wires, best_payload

    def connections(self, parts):
        last, pairs = {}, defaultdict(list)
        for i in parts:
            atom = self.occurrences[i]
            for port, q in enumerate(atom.qubits):
                if q in last:
                    previous, previous_port = last[q]
                    pairs[previous, i].append((previous_port, port))
                last[q] = i, port
        return [(a, b, tuple(sorted(ports))) for (a, b), ports in sorted(pairs.items())]

    def register(self, parts, round_number, parent=None, added=None):
        self.check()
        gates = frozenset(g for i in parts for g in self.occurrences[i].gates)
        sig, wires, payload = self.signature(parts)
        key = gates, sig
        if key in self.region_keys:
            return self.region_keys[key], False
        if len(self.regions) >= self.limits.max_regions:
            raise _Stopped("regions_limit")
        if sig in self.signatures and self.signatures[sig] != payload:
            raise RuntimeError("signature collision: unequal structures have the same digest")
        pending = []
        for a, b, ports in self.connections(parts):
            left, right = self.occurrences[a], self.occurrences[b]
            link_key = left.class_id, right.class_id, ports
            witness = left.id, right.id
            existing = self.link_keys.get(link_key)
            if existing is None or witness not in self.links[existing]["_witnesses"]:
                pending.append((link_key, witness))
        if self.counts["link_witnesses"] + len(pending) > self.limits.max_link_witnesses:
            raise _Stopped("link_witnesses_limit")
        for link_key, witness in pending:
            if link_key not in self.link_keys:
                self.link_keys[link_key] = len(self.links)
                self.links.append({"id": f"l{len(self.links)}", "source_class": link_key[0],
                                   "target_class": link_key[1], "ports": link_key[2],
                                   "_witnesses": set()})
            self.links[self.link_keys[link_key]]["_witnesses"].add(witness)
            self.counts["link_witnesses"] += 1
        prefix, suffix = self.dag.boundary(gates)
        region = _Region(f"r{len(self.regions)}", parts, gates, wires, sig,
                         prefix, suffix, round_number, parent, added)
        self.region_keys[key] = len(self.regions)
        self.signatures[sig] = payload
        self.regions.append(region)
        return len(self.regions) - 1, True

    def candidates(self, region, side):
        first, last = self.ports(region.gates)
        found = set()
        if side == "prefix":
            for q, gid in last.items():
                nxt = self.dag.next.get((gid, q))
                found.update(self.starts.get((nxt, q), ()))
        else:
            for q, gid in first.items():
                prev = self.dag.prev.get((gid, q))
                found.update(self.ends.get((prev, q), ()))
        return sorted(found)

    def compatible(self, region, occurrence, side):
        if region.gates & occurrence.gates:
            return False
        union = region.gates | occurrence.gates
        if len(union) > self.limits.max_gates:
            self.counts["size_pruned"] += 1
            return False
        rf, rl = self.ports(region.gates)
        af, al = self.ports(occurrence.gates)
        if len(rf.keys() | af.keys()) > self.limits.max_qubits:
            self.counts["width_pruned"] += 1
            return False
        shared = rf.keys() & af.keys()
        if not shared:
            return False
        for q in shared:
            left, right = (rl[q], af[q]) if side == "prefix" else (al[q], rf[q])
            if self.dag.next.get((left, q)) != right:
                return False
        # Requiring predecessor/successor closure *at each intermediate step*
        # loses valid fan-in/out regions. Require convexity here; compute the
        # actual prefix/suffix qualifications separately for the completed union.
        return self.convex(union)

    def snapshot(self, number, stage, before_regions, before_links):
        info = {"round": number, "stage": stage, "matches": len(self.occurrences),
                "regions": len(self.regions), "new_regions": len(self.regions) - before_regions,
                "links": len(self.links), "new_links": len(self.links) - before_links,
                "elapsed_seconds": time.perf_counter() - self.start}
        self.rounds.append(info)
        if self.on_round:
            self.on_round(dict(info))

    def run(self):
        reason = "fixed_point_within_bounds"
        stage, round_number, before_regions, before_links = "matching", 0, 0, 0
        try:
            self.seed()
            self.snapshot(0, stage, 0, 0)
            frontier = [(i, side) for i in range(len(self.regions))
                        for side in ("prefix", "suffix")]
            stage = "linking"
            while frontier:
                self.check()
                if round_number >= self.limits.max_rounds:
                    raise _Stopped("rounds_limit")
                round_number += 1
                before_regions, before_links = len(self.regions), len(self.links)
                next_frontier = []
                for region_id, side in frontier:
                    region = self.regions[region_id]
                    for occurrence_id in self.candidates(region, side):
                        self.check("extension_attempts")
                        occurrence = self.occurrences[occurrence_id]
                        if not self.compatible(region, occurrence, side):
                            continue
                        parts = (region.parts + (occurrence_id,)) if side == "prefix" else (
                            (occurrence_id,) + region.parts)
                        new_id, created = self.register(parts, round_number, region.id, occurrence.id)
                        if created:
                            for direction in ("prefix", "suffix"):
                                next_frontier.append((new_id, direction))
                self.snapshot(round_number, stage, before_regions, before_links)
                frontier = next_frontier
        except _Stopped as stopped:
            reason = str(stopped)
            # Include partial progress without duplicating a completed round.
            if not self.rounds or self.rounds[-1]["round"] != round_number or stage == "matching":
                self.snapshot(round_number, stage, before_regions, before_links)
        return self.result(reason)

    def result(self, reason):
        matches, regions, groups = [], [], defaultdict(list)
        annotations = {gid: {"matches": [], "regions": []} for gid in self.dag.order}
        for o in self.occurrences:
            gates = sorted(o.gates, key=self.dag.rank.__getitem__)
            matches.append({"id": o.id, "class_id": o.class_id, "member_id": o.member_id,
                            "gate_ids": gates, "qubit_map": list(o.qubits), "bindings": dict(o.bindings),
                            "boundary": {"prefix": o.prefix, "suffix": o.suffix}})
            for g in gates:
                annotations[g]["matches"].append(o.id)
        for r in self.regions:
            groups[r.signature].append(r.id)
            first, last = self.ports(r.gates)
            regions.append({"id": r.id, "composition_key": r.signature,
                            "parts": [self.occurrences[i].id for i in r.parts],
                            "gate_ids": sorted(r.gates, key=self.dag.rank.__getitem__),
                            "qubit_map": list(r.qubits), "round": r.round,
                            "parent_region": r.parent, "added_match": r.added,
                            "boundary": {"prefix": r.prefix, "suffix": r.suffix},
                            "input_cut": [{"qubit": q, "first_gate": first[q],
                                           "predecessor": self.dag.prev.get((first[q], q))}
                                          for q in sorted(first)],
                            "output_cut": [{"qubit": q, "last_gate": last[q],
                                            "successor": self.dag.next.get((last[q], q))}
                                           for q in sorted(last)]})
            for g in r.gates:
                annotations[g]["regions"].append(r.id)
        return {"schema": "qclassification/1",
                "dag_fingerprint": digest(self.dag.to_dict()),
                "library_fingerprint": digest(self.library.to_dict()),
                "global_phase": self.dag.global_phase,
                "catalog_invariant": {"classes_before": len(self.library.classes),
                                      "classes_after": len(self.library.classes),
                                      "members_before": sum(len(c.members) for c in self.library.classes),
                                      "members_after": sum(len(c.members) for c in self.library.classes)},
                "stop_reason": reason, "complete_for_bounds": reason == "fixed_point_within_bounds",
                "limits": asdict(self.limits), "rounds": self.rounds,
                "statistics": {**self.counts, "matches": len(matches), "regions": len(regions),
                               "links": len(self.links), "composition_groups": len(groups),
                               "prefix_regions": sum(r.prefix for r in self.regions),
                               "suffix_regions": sum(r.suffix for r in self.regions),
                               "elapsed_seconds": time.perf_counter() - self.start},
                "matches": matches,
                "links": [{"id": link["id"], "source_class": link["source_class"],
                           "target_class": link["target_class"], "ports": [list(p) for p in link["ports"]],
                           "witnesses": [{"source_match": a, "target_match": b}
                                         for a, b in sorted(link["_witnesses"])]} for link in self.links],
                "regions": regions,
                "composition_groups": [{"key": sig, "signature": self.signatures[sig],
                                        "region_ids": ids} for sig, ids in sorted(groups.items())],
                "annotations": annotations,
                "unmatched_gate_ids": [g for g in self.dag.order if not annotations[g]["matches"]],
                "notes": ["Class/member definitions and circuit are unchanged.",
                          "Composition keys describe link structures, not new equivalence classes.",
                          "Groups may overlap; no reuse selection or extraction has been performed.",
                          "Suffix and interior operator matches do not imply zero-input state reuse."]}


def classify(dag: CircuitDAG, library: ClassLibrary, limits=None, *, on_round=None,
             enable_overlap=False, overlap_limits=None):
    """Classify a supplied DAG using a trusted, fixed supplied class library.

    Return a JSON-compatible result. Budget exhaustion returns a sound partial
    result, distinguished from completion by stop_reason/complete_for_bounds.
    Callbacks receive one progress dictionary per completed or interrupted round.
    With enable_overlap=True, add independent bounded-search status in 'overlap'
    and whole-window alternatives in 'rewrite_candidates'. Original labels and
    top-level completion status continue to describe the baseline classifier.
    """
    if type(enable_overlap) is not bool:
        raise InputError("enable_overlap must be a boolean")
    if overlap_limits is not None and not enable_overlap:
        raise InputError("overlap_limits requires enable_overlap=True")
    if enable_overlap:
        from .overlap import OverlapLimits, search_overlaps
        overlap_limits = overlap_limits or OverlapLimits()
        overlap_limits.validate()
    limits = limits or Limits()
    limits.validate()
    result = _Run(dag, library, limits, on_round).run()
    if library.equivalence == "up_to_global_phase":
        result["equivalence"] = library.equivalence
        result["notes"].append(
            "Library equality is exact modulo global phase. Relative global phases are not tracked; "
            "grouped operators/states need not have identical complex amplitudes.")
    if enable_overlap:
        result["overlap"], result["rewrite_candidates"] = search_overlaps(dag, library, result, overlap_limits)
    return result
