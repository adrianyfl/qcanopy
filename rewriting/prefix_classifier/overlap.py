"""Bounded two-step, partial-RHS overlap search. No input graph is mutated."""
from collections import defaultdict, deque
from dataclasses import asdict, dataclass
import heapq
import time

from .engine import Limits, _Run
from .model import CircuitDAG, InputError, Variable, canonical_json, digest, param_json


@dataclass(frozen=True)
class OverlapLimits:
    max_qubits: int = 4
    max_original_gates: int = 12
    max_virtual_gates: int = 16
    max_context_steps: int = 3
    max_windows_per_match: int = 8
    max_members_per_class: int = 4
    max_views_per_match: int = 8
    max_views: int = 512
    max_candidates: int = 256
    max_match_attempts: int = 50000
    max_graph_steps: int = 100000
    time_limit_seconds: float = 2.0
    seed_order: str = "boundary_distance"
    focus: str = "both"
    max_candidates_per_window: int = 64

    def validate(self):
        for name, value in asdict(self).items():
            if name == "seed_order":
                if value not in ("boundary_distance", "library"):
                    raise InputError("overlap seed_order must be boundary_distance or library")
            elif name == "focus":
                if value not in ("prefix", "suffix", "both"):
                    raise InputError("overlap focus must be prefix, suffix or both")
            elif name == "time_limit_seconds":
                if type(value) not in (int, float) or not (0 < value < float("inf")):
                    raise InputError("overlap time_limit_seconds must be positive and finite")
            elif type(value) is not int or value < 1:
                raise InputError(f"overlap {name}: expected a positive integer")
        if self.max_qubits > 6 or max(self.max_original_gates, self.max_virtual_gates) > 256:
            raise InputError("overlap caps width at 6 and window sizes at 256")


class _Stop(Exception):
    pass


class _WindowFull(Exception):
    pass


def _shape(gate):
    return gate.op, len(gate.qubits), len(gate.params)


def _convex_hull(dag, gates, check):
    """All vertices on paths between selected gates, within their rank interval."""
    low = min(dag.rank[g] for g in gates)
    high = max(dag.rank[g] for g in gates)

    def walk(adjacency):
        seen, todo = set(), list(gates)
        while todo:
            check("graph_steps")
            g = todo.pop()
            if g in seen or not low <= dag.rank[g] <= high:
                continue
            seen.add(g)
            todo.extend(adjacency[g])
        return seen

    return frozenset(walk(dag.successors) & walk(dag.predecessors))


def _window(dag, gates):
    """Extract a convex window with local port numbering and zero local phase."""
    wires = sorted({q for g in gates for q in dag.by_id[g].qubits})
    mapping = {q: i for i, q in enumerate(wires)}
    data = {"schema": "qdag/1", "num_qubits": len(wires), "global_phase": 0,
            "nodes": [{"id": g, "op": dag.by_id[g].op,
                       "qubits": [mapping[q] for q in dag.by_id[g].qubits],
                       "params": [param_json(p) for p in dag.by_id[g].params]}
                      for g in sorted(gates)],
            "edges": [{"source": dag.prev[b, q], "target": b, "qubit": mapping[q]}
                      for b in sorted(gates) for q in dag.by_id[b].qubits
                      if dag.prev.get((b, q)) in gates]}
    return CircuitDAG.from_dict(data), wires


def _replace(dag, consumed, member, qubits, bindings, prefix):
    """Splice a whole same-interface member at a convex occurrence's wire cuts."""
    consumed = frozenset(consumed)
    first, last = dag.ports(consumed)
    # Even adversarial input IDs cannot alias the view's virtual IDs.
    while any(g.startswith(prefix) for g in dag.by_id):
        prefix += "_"
    ids = {g: prefix + g for g in member.dag.order}
    replacement = member.dag.to_dict()
    for node in replacement["nodes"]:
        gate = member.dag.by_id[node["id"]]
        node["id"] = ids[gate.id]
        node["qubits"] = [qubits[q] for q in gate.qubits]
        node["params"] = [param_json(bindings[p.name] if isinstance(p, Variable) else p)
                          for p in gate.params]
    data = dag.to_dict()
    data["nodes"] = [g for g in data["nodes"] if g["id"] not in consumed] + replacement["nodes"]
    data["edges"] = [e for e in data["edges"]
                     if e["source"] not in consumed and e["target"] not in consumed]
    data["edges"].extend({"source": ids[a], "target": ids[b], "qubit": qubits[q]}
                         for a, b, q in member.dag.edges)
    mf, ml = member.dag.ports(frozenset(member.dag.by_id))
    for port, q in enumerate(qubits):
        prev, nxt = dag.prev.get((first[q], q)), dag.next.get((last[q], q))
        if prev is not None:
            data["edges"].append({"source": prev, "target": ids[mf[port]], "qubit": q})
        if nxt is not None:
            data["edges"].append({"source": ids[ml[port]], "target": nxt, "qubit": q})
    return CircuitDAG.from_dict(data), tuple(ids[g] for g in member.dag.order)


def _structure(dag):
    """ID-independent structure, retaining ordered local ports and exact parameters."""
    labels = {g.id: {"op": g.op, "qubits": list(g.qubits),
                     "params": [param_json(p) for p in g.params]} for g in dag.nodes}
    degrees = {g: len(dag.predecessors[g]) for g in dag.order}
    ready = [(canonical_json(labels[g]), g) for g in dag.order if not degrees[g]]
    heapq.heapify(ready)
    ordered = []
    while ready:
        _, g = heapq.heappop(ready)
        ordered.append(labels[g])
        for nxt in dag.successors[g]:
            degrees[nxt] -= 1
            if not degrees[nxt]:
                heapq.heappush(ready, (canonical_json(labels[nxt]), nxt))
    return canonical_json({"num_qubits": dag.num_qubits, "gates": ordered})


class _ViewMatcher(_Run):
    """Reuse the exact matcher, charging every view to the same overlap budget."""
    def __init__(self, dag, owner):
        self.owner = owner
        super().__init__(dag, owner.library,
                         Limits(max_qubits=owner.limits.max_qubits,
                                max_gates=owner.limits.max_virtual_gates), None)

    def check(self, counter=None):
        self.owner.check(counter)


class _OverlapSearch:
    def __init__(self, dag, library, base, limits):
        self.dag, self.library, self.base, self.limits = dag, library, base, limits
        self.start = time.perf_counter()
        self.counts = {name: 0 for name in (
            "views", "match_attempts", "graph_steps", "windows", "partial_overlaps",
            "window_pruned", "member_pruned", "view_pruned", "target_members_pruned",
            "duplicate_candidates", "unchanged_candidates", "candidate_members_examined",
            "window_candidate_pruned", "seeds_visited")}
        self.by_class = {c.id: c for c in library.classes}
        self.member_index = defaultdict(set)
        self.members = {}
        self.candidates, self.keys = [], set()
        self.window_counts = defaultdict(int)

    def check(self, counter=None):
        if time.perf_counter() - self.start >= self.limits.time_limit_seconds:
            raise _Stop("time_limit")
        if counter:
            if self.counts[counter] >= getattr(self.limits, "max_" + counter):
                raise _Stop(counter + "_limit")
            self.counts[counter] += 1

    def index(self):
        for c in self.library.classes:
            for m in c.members:
                self.check()
                if c.num_qubits > self.limits.max_qubits or len(m.dag.nodes) > self.limits.max_virtual_gates:
                    self.counts["target_members_pruned"] += 1
                    continue
                key = c.id, m.id
                self.members[key] = c, m
                for g in m.dag.nodes:
                    self.check()
                    self.member_index[_shape(g)].add(key)

    def fits(self, gates):
        if len(gates) > self.limits.max_original_gates:
            return False
        return len({q for g in gates for q in self.dag.by_id[g].qubits}) <= self.limits.max_qubits

    def windows(self, seed):
        """BFS over convex windows formed by adding one adjacent gate then closing paths."""
        if not self.fits(seed):
            self.counts["window_pruned"] += 1
            return
        seen, todo, emitted = {seed}, deque([(seed, 0)]), 0
        while todo:
            self.check()
            current, depth = todo.popleft()
            if depth >= self.limits.max_context_steps:
                continue
            neighbors = {h for g in current for h in
                         self.dag.predecessors[g] + self.dag.successors[g] if h not in current}
            for neighbor in sorted(neighbors, key=self.dag.rank.__getitem__):
                self.check()
                enlarged = _convex_hull(self.dag, current | {neighbor}, self.check)
                if enlarged in seen:
                    continue
                if not self.fits(enlarged):
                    self.counts["window_pruned"] += 1
                    continue
                if emitted >= self.limits.max_windows_per_match:
                    self.counts["window_pruned"] += 1
                    return
                seen.add(enlarged)
                emitted += 1
                self.counts["windows"] += 1
                yield enlarged
                todo.append((enlarged, depth + 1))

    def alternatives(self, eqclass, matched_member, *, partial=False):
        members = [m for m in eqclass.members if m.id != matched_member and
                   (not partial or len(m.dag.nodes) >= 2)]
        if len(members) > self.limits.max_members_per_class:
            self.counts["member_pruned"] += len(members) - self.limits.max_members_per_class
        return members[:self.limits.max_members_per_class]

    def target_members(self, view, introduced):
        # A crossing match must contain at least one virtual and one untouched gate.
        inside, outside = set(), set()
        for g in view.nodes:
            self.check()
            (inside if g.id in introduced else outside).update(self.member_index[_shape(g)])
        return [self.members[k] for k in sorted(inside & outside)]

    def save(self, occurrence, source_member, original, wires, support,
             intermediate, introduced, target_class, target_member, matched,
             target_wires, bindings, replacement, final, added):
        structure = _structure(final)
        if structure == _structure(original):
            self.counts["unchanged_candidates"] += 1
            return
        key = support, structure
        if key in self.keys:
            self.counts["duplicate_candidates"] += 1
            return
        if self.window_counts[support] >= self.limits.max_candidates_per_window:
            self.counts["window_candidate_pruned"] += 1
            raise _WindowFull()
        if len(self.candidates) >= self.limits.max_candidates:
            raise _Stop("candidates_limit")
        prefix, suffix = self.dag.boundary(support)
        first, last = self.dag.ports(support)
        consumed = sorted(matched, key=intermediate.rank.__getitem__)
        self.candidates.append({
            "id": f"rw{len(self.candidates)}", "origin_match_id": occurrence["id"],
            "equivalence": self.library.equivalence,
            "original_gate_ids": sorted(support, key=self.dag.rank.__getitem__),
            "qubit_map": wires,
            "boundary": {"prefix": prefix, "suffix": suffix},
            "input_cut": [{"qubit": q, "first_gate": first[q],
                           "predecessor": self.dag.prev.get((first[q], q))} for q in wires],
            "output_cut": [{"qubit": q, "last_gate": last[q],
                            "successor": self.dag.next.get((last[q], q))} for q in wires],
            "original_dag": original.to_dict(),
            "intermediate_dag": intermediate.to_dict(), "alternative_dag": final.to_dict(),
            "alternative_key": digest(structure),
            "contextual_match": {
                "view": "intermediate_dag", "class_id": target_class.id,
                "member_id": target_member.id, "gate_ids": consumed,
                "qubit_map": list(target_wires), "bindings": bindings,
                "consumed_virtual_gate_ids": [g for g in consumed if g in introduced],
                "original_context_gate_ids": [g for g in consumed if g not in introduced],
                "retained_virtual_gate_ids": [g for g in introduced if g not in matched]},
            "proof": [
                {"class_id": occurrence["class_id"], "from_member_id": occurrence["member_id"],
                 "to_member_id": source_member.id, "consumed_gate_ids": occurrence["gate_ids"],
                 "introduced_gate_ids": list(introduced),
                 "qubit_map": [wires.index(q) for q in occurrence["qubit_map"]],
                 "bindings": occurrence["bindings"]},
                {"class_id": target_class.id, "from_member_id": target_member.id,
                 "to_member_id": replacement.id, "consumed_gate_ids": consumed,
                 "introduced_gate_ids": list(added), "qubit_map": list(target_wires),
                 "bindings": bindings}],
        })
        self.keys.add(key)
        self.window_counts[support] += 1

    def search_occurrence(self, occurrence):
        eqclass = self.by_class[occurrence["class_id"]]
        alternatives = self.alternatives(eqclass, occurrence["member_id"], partial=True)
        if not alternatives:
            return
        seed = frozenset(occurrence["gate_ids"])
        views = 0
        for support in self.windows(seed):
            if self.window_counts[support] >= self.limits.max_candidates_per_window:
                self.counts["window_candidate_pruned"] += 1
                continue
            before = self.counts["views"]
            try:
                self.search_window(occurrence, support, alternatives,
                                   self.limits.max_views_per_match - views)
            except _WindowFull:
                pass
            views += self.counts["views"] - before
            if views >= self.limits.max_views_per_match:
                self.counts["view_pruned"] += 1
                return

    def search_window(self, occurrence, support, alternatives, view_budget):
        seed = frozenset(occurrence["gate_ids"])
        original, wires = _window(self.dag, support)
        source_wires = tuple(wires.index(q) for q in occurrence["qubit_map"])
        views = 0
        for member in alternatives:
            self.check()
            if len(support) - len(seed) + len(member.dag.nodes) > self.limits.max_virtual_gates:
                self.counts["view_pruned"] += 1
                continue
            if views >= view_budget:
                self.counts["view_pruned"] += 1
                return
            self.check("views")
            views += 1
            view, introduced = _replace(original, seed, member, source_wires,
                                        occurrence["bindings"], "@overlap1:")
            introduced_set = frozenset(introduced)
            matcher = _ViewMatcher(view, self)
            for target_class, target_member in self.target_members(view, introduced_set):
                self.check()
                self.counts["candidate_members_examined"] += 1
                for matched, target_wires, binding_pairs in matcher.match_member(target_class, target_member):
                    intersection = matched & introduced_set
                    if not intersection or intersection == introduced_set or not (matched - introduced_set):
                        continue
                    self.counts["partial_overlaps"] += 1
                    bindings = dict(binding_pairs)
                    for replacement in self.alternatives(target_class, target_member.id):
                        self.check()
                        if len(view.nodes) - len(matched) + len(replacement.dag.nodes) > self.limits.max_virtual_gates:
                            self.counts["view_pruned"] += 1
                            continue
                        final, added = _replace(view, matched, replacement, target_wires,
                                                bindings, "@overlap2:")
                        self.save(occurrence, member, original, wires, support, view,
                                  introduced, target_class, target_member, matched,
                                  target_wires, bindings, replacement, final, added)

    def scheduled_occurrences(self):
        if self.limits.seed_order == "library":
            return self.base["matches"]
        # Distances use only the graph. No truth tags, target IDs or rule IDs.
        starts, ends = {}, {}
        wires = defaultdict(list)
        for gid in self.dag.order:
            for q in self.dag.by_id[gid].qubits:
                wires[q].append(gid)
        for ids in wires.values():
            for i, gid in enumerate(ids):
                starts[gid] = min(starts.get(gid, i), i)
                back = len(ids) - i - 1
                ends[gid] = min(ends.get(gid, back), back)
        def ordered(distance):
            return sorted(self.base["matches"],
                          key=lambda m: min(distance[g] for g in m["gate_ids"]))
        if self.limits.focus != "both":
            return ordered(starts if self.limits.focus == "prefix" else ends)
        # Alternate the two boundary queues so one side cannot take all slots.
        front, back = ordered(starts), ordered(ends)
        out, seen = [], set()
        for a, b in zip(front, back):
            for m in (a, b):
                if m["id"] not in seen:
                    out.append(m); seen.add(m["id"])
        return out

    def run(self):
        reason = "completed_bounded_search"
        try:
            self.index()
            for occurrence in self.scheduled_occurrences():
                self.check()
                self.counts["seeds_visited"] += 1
                self.search_occurrence(occurrence)
        except _Stop as stopped:
            reason = str(stopped)
        if reason == "completed_bounded_search":
            if any(self.counts[k] for k in ("window_pruned", "member_pruned", "view_pruned", "target_members_pruned", "window_candidate_pruned")):
                reason = "search_pruned"
            elif not self.base["complete_for_bounds"]:
                reason = "base_classification_incomplete"
        return {
            "enabled": True, "stop_reason": reason,
            "equivalence": self.library.equivalence,
            "complete_for_bounds": reason == "completed_bounded_search",
            "base_complete_for_bounds": self.base["complete_for_bounds"],
            "scope": "one member substitution followed by one proper partial-RHS crossing substitution",
            "limits": asdict(self.limits),
            "statistics": {**self.counts, "candidates": len(self.candidates),
                           "elapsed_seconds": time.perf_counter() - self.start},
            "notes": ["Candidates are alternative whole-window plans, not original-DAG class labels.",
                      "Local DAGs have zero explicit phase; retain the full input circuit global_phase.",
                      "Candidate boundary flags describe the whole original window, not an internal virtual match.",
                      "For up_to_global_phase libraries, relative rewrite phases are not tracked.",
                      "No extraction, recursive expansion, class union or cross-view composition is performed."]}, self.candidates


def search_overlaps(dag, library, base, limits):
    return _OverlapSearch(dag, library, base, limits).run()
