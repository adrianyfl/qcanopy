import cmath
import copy
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from prefix_classifier import CircuitDAG, ClassLibrary, InputError, Limits, OverlapLimits, classify


ROOT = Path(__file__).resolve().parents[1]


def dag_data(n, steps, phase=0):
    nodes, edges, last = [], [], {}
    for i, (op, qs, *params) in enumerate(steps):
        gid = f"g{i}"
        nodes.append({"id": gid, "op": op, "qubits": qs, "params": params[0] if params else []})
        for q in qs:
            if q in last:
                edges.append({"source": last[q], "target": gid, "qubit": q})
            last[q] = gid
    return {"schema": "qdag/1", "num_qubits": n, "nodes": nodes,
            "edges": edges, "global_phase": phase}


def eqclass(cid, n, members, parameters=()):
    return {"id": cid, "num_qubits": n, "parameters": list(parameters),
            "members": [{"id": name, "dag": dag_data(n, steps)} for name, steps in members]}


def library(classes=None):
    if classes is None:
        classes = [eqclass("commute", 1, [
            ("zh", [("z", [0]), ("h", [0])]),
            ("hx", [("h", [0]), ("x", [0])])]),
            eqclass("cancel_x", 1, [("xx", [("x", [0]), ("x", [0])]),
                                     ("id", [("id", [0])])])]
    return ClassLibrary.from_dict({"schema": "qclass-library/1", "equivalence": "exact_unitary",
                                   "classes": classes})


def simulate(dag, initial):
    """Independent amplitude calculation on the DAG, not on rewrite proofs."""
    state = list(initial)
    for gid in dag.order:
        gate = dag.by_id[gid]
        if gate.op == "id":
            continue
        out = [0j] * len(state)
        if gate.op == "cx":
            c, t = gate.qubits
            for i, value in enumerate(state):
                j = i ^ (1 << t) if (i >> c) & 1 else i
                out[j] += value
        else:
            q = gate.qubits[0]
            if gate.op == "h":
                s = 1 / math.sqrt(2)
                matrix = ((s, s), (s, -s))
            elif gate.op == "x":
                matrix = ((0, 1), (1, 0))
            elif gate.op == "z":
                matrix = ((1, 0), (0, -1))
            elif gate.op == "rz":
                theta = gate.params[0]
                matrix = ((cmath.exp(-0.5j * theta), 0), (0, cmath.exp(0.5j * theta)))
            else:
                raise AssertionError(gate.op)
            for i, value in enumerate(state):
                for row in (0, 1):
                    out[(i & ~(1 << q)) | (row << q)] += matrix[row][(i >> q) & 1] * value
        state = out
    return [cmath.exp(1j * dag.global_phase) * value for value in state]


def stable_base(result):
    result = copy.deepcopy(result)
    result.pop("overlap", None)
    result.pop("rewrite_candidates", None)
    result["statistics"].pop("elapsed_seconds")
    for row in result["rounds"]:
        row.pop("elapsed_seconds")
    return result


class OverlapTests(unittest.TestCase):
    def test_boundary_scheduling_and_independent_side_focus(self):
        from prefix_classifier.overlap import _OverlapSearch
        d = CircuitDAG.from_dict(dag_data(1, [("h", [0]) for _ in range(6)]))
        ms = [dict(id=name, gate_ids=[gid]) for name,gid in
              [('middle','g2'),('tail','g5'),('head','g0')]]
        for focus,expected in [('prefix','head'),('suffix','tail'),('both','head')]:
            search = _OverlapSearch(d,library(),{'matches':ms},OverlapLimits(focus=focus))
            actual = search.scheduled_occurrences()
            self.assertEqual(actual[0]['id'],expected)
            self.assertEqual({m['id'] for m in actual},{'middle','tail','head'})
            if focus=='both':self.assertEqual(actual[1]['id'],'tail')
        legacy = _OverlapSearch(d,library(),{'matches':ms},OverlapLimits(seed_order='library'))
        self.assertEqual(legacy.scheduled_occurrences(),ms)

    def test_physical_window_quota_limits_alternatives_not_soundness(self):
        lib=library([eqclass('commute',1,[('zh',[('z',[0]),('h',[0])]),
                                        ('hx',[('h',[0]),('x',[0])])]),
                     eqclass('cancel_x',1,[('xx',[('x',[0]),('x',[0])]),
                         ('id',[('id',[0])]),('hh',[('h',[0]),('h',[0])]),
                         ('zz',[('z',[0]),('z',[0])])])])
        r=classify(self.dag,lib,enable_overlap=True,
                   overlap_limits=OverlapLimits(max_candidates_per_window=1))
        from collections import Counter
        sizes=Counter(tuple(c['original_gate_ids']) for c in r['rewrite_candidates'])
        self.assertTrue(sizes)
        self.assertLessEqual(max(sizes.values()),1)
        self.assertGreater(r['overlap']['statistics']['window_candidate_pruned'],0)
        self.assertFalse(r['overlap']['complete_for_bounds'])
        for c in r['rewrite_candidates']:self.assert_sound(c)

    def setUp(self):
        self.dag = CircuitDAG.from_dict(dag_data(1, [("z", [0]), ("h", [0]), ("x", [0])], 0.37))
        self.lib = library()

    def assert_sound(self, candidate):
        dags = [CircuitDAG.from_dict(candidate[k])
                for k in ("original_dag", "intermediate_dag", "alternative_dag")]
        n = dags[0].num_qubits
        for basis in range(1 << n):
            state = [complex(i == basis) for i in range(1 << n)]
            expected = simulate(dags[0], state)
            for d in dags[1:]:
                actual = simulate(d, state)
                self.assertLess(max(abs(a - b) for a, b in zip(actual, expected)), 1e-12)
        # Proof step IDs refer to the preceding view, never to a guessed B->Y mapping.
        for i, proof in enumerate(candidate["proof"]):
            consumed, added = set(proof["consumed_gate_ids"]), set(proof["introduced_gate_ids"])
            self.assertTrue(consumed <= set(dags[i].by_id))
            self.assertTrue(added <= set(dags[i + 1].by_id))
            self.assertEqual(set(dags[i + 1].by_id), (set(dags[i].by_id) - consumed) | added)

    def test_disabled_is_legacy_and_explicit_flag_is_required(self):
        with patch("prefix_classifier.overlap.search_overlaps", side_effect=AssertionError("called")):
            result = classify(self.dag, self.lib)
            self.assertNotIn("overlap", result)
            self.assertNotIn("rewrite_candidates", result)
        with self.assertRaises(InputError):
            classify(self.dag, self.lib, overlap_limits=OverlapLimits())
        with self.assertRaises(InputError):
            classify(self.dag, self.lib, enable_overlap="yes")

    def test_partial_rhs_keeps_residual_and_does_not_pollute_original_matches(self):
        before = self.dag.to_dict(), self.lib.to_dict()
        base = classify(self.dag, self.lib)
        result = classify(self.dag, self.lib, enable_overlap=True)
        self.assertEqual(stable_base(base), stable_base(result))
        self.assertEqual(before, (self.dag.to_dict(), self.lib.to_dict()))
        self.assertTrue(result["overlap"]["complete_for_bounds"])
        self.assertEqual(len(result["rewrite_candidates"]), 1)
        c = result["rewrite_candidates"][0]
        self.assertEqual(c["original_gate_ids"], ["g0", "g1", "g2"])
        self.assertEqual(c["proof"][0]["from_member_id"], "zh")
        self.assertEqual(c["proof"][0]["to_member_id"], "hx")
        self.assertEqual(c["proof"][1]["from_member_id"], "xx")
        self.assertEqual(c["proof"][1]["to_member_id"], "id")
        self.assertEqual(c["contextual_match"]["original_context_gate_ids"], ["g2"])
        self.assertEqual(len(c["contextual_match"]["retained_virtual_gate_ids"]), 1)
        self.assertFalse(any(m["class_id"] == "cancel_x" for m in result["matches"]))
        final = CircuitDAG.from_dict(c["alternative_dag"])
        self.assertEqual([final.by_id[g].op for g in final.order], ["h", "id"])
        self.assertEqual(result["global_phase"], 0.37)
        self.assertEqual(c["alternative_dag"]["global_phase"], 0)
        self.assert_sound(c)

    def test_multiwire_entangled_case_and_nontrivial_port_map(self):
        lib = library([eqclass("move_x", 2, [
            ("before", [("x", [0]), ("cx", [0, 1])]),
            ("after", [("cx", [0, 1]), ("x", [0]), ("x", [1])])]),
            eqclass("cancel_x", 1, [("xx", [("x", [0]), ("x", [0])]), ("id", [("id", [0])])])])
        dag = CircuitDAG.from_dict(dag_data(3, [("x", [2]), ("cx", [2, 0]), ("x", [0])]))
        result = classify(dag, lib, enable_overlap=True)
        self.assertTrue(result["rewrite_candidates"])
        for c in result["rewrite_candidates"]:
            self.assertEqual(c["qubit_map"], [0, 2])
            self.assertEqual(c["proof"][0]["qubit_map"], [1, 0])
            self.assertEqual(len(c["contextual_match"]["retained_virtual_gate_ids"]), 2)
            self.assert_sound(c)

    def test_no_false_match_with_reversed_control_target(self):
        lib = library([eqclass("commute_z", 2, [
            ("before", [("z", [0]), ("cx", [0, 1])]),
            ("after", [("cx", [0, 1]), ("z", [0])])])])
        dag = CircuitDAG.from_dict(dag_data(2, [("z", [1]), ("cx", [1, 0]), ("cx", [0, 1])]))
        result = classify(dag, lib, enable_overlap=True)
        self.assertFalse(result["rewrite_candidates"])

    def test_incoming_boundary_overlap_keeps_right_residual(self):
        lib = library([eqclass("move", 1, [
            ("hz", [("h", [0]), ("z", [0])]),
            ("xh", [("x", [0]), ("h", [0])])]),
            eqclass("cancel_x", 1, [("xx", [("x", [0]), ("x", [0])]), ("id", [("id", [0])])])])
        dag = CircuitDAG.from_dict(dag_data(1, [("x", [0]), ("h", [0]), ("z", [0])]))
        result = classify(dag, lib, enable_overlap=True)
        self.assertTrue(result["rewrite_candidates"])
        c = result["rewrite_candidates"][0]
        self.assertEqual(c["contextual_match"]["original_context_gate_ids"], ["g0"])
        final = CircuitDAG.from_dict(c["alternative_dag"])
        self.assertEqual([final.by_id[g].op for g in final.order], ["id", "h"])
        self.assert_sound(c)

    def test_whole_replacement_is_not_reported_as_partial_overlap(self):
        lib = library([eqclass("commute", 1, [
            ("zh", [("z", [0]), ("h", [0])]), ("hx", [("h", [0]), ("x", [0])])]),
            eqclass("whole", 1, [("hxx", [("h", [0]), ("x", [0]), ("x", [0])]),
                                    ("h", [("h", [0])])])])
        result = classify(self.dag, lib, enable_overlap=True)
        # The reverse h->hxx direction can legitimately create OTHER partial
        # overlaps. Only the whole hx + outside x -> h case must be absent.
        self.assertFalse(any(c["proof"][0]["class_id"] == "commute" and
                             c["proof"][1]["class_id"] == "whole"
                             for c in result["rewrite_candidates"]))
        for c in result["rewrite_candidates"]:
            self.assert_sound(c)

    def test_parameters_bind_and_conflicts_are_rejected(self):
        t = {"var": "theta"}
        lib = library([eqclass("move_rz", 1, [
            ("before", [("rz", [0], [t]), ("z", [0])]),
            ("after", [("z", [0]), ("rz", [0], [t])])], ["theta"]),
            eqclass("insert_zz", 1, [
                ("pair", [("rz", [0], [t]), ("rz", [0], [t])]),
                ("with_zz", [("rz", [0], [t]), ("z", [0]), ("z", [0]), ("rz", [0], [t])])], ["theta"])])
        for theta, expected in ((0.25, True), (0.5, False)):
            dag = CircuitDAG.from_dict(dag_data(1, [("rz", [0], [0.25]), ("z", [0]), ("rz", [0], [theta])]))
            result = classify(dag, lib, enable_overlap=True)
            candidates = [c for c in result["rewrite_candidates"] if c["proof"][1]["class_id"] == "insert_zz"]
            self.assertEqual(bool(candidates), expected)
            for c in candidates:
                self.assertEqual(c["proof"][1]["bindings"], {"theta": 0.25})
                self.assert_sound(c)

    def test_input_array_order_and_virtual_id_collisions(self):
        data = self.dag.to_dict()
        renames = {g["id"]: "@overlap1:" + g["id"] for g in data["nodes"]}
        for node in data["nodes"]:
            node["id"] = renames[node["id"]]
        for edge in data["edges"]:
            edge["source"], edge["target"] = renames[edge["source"]], renames[edge["target"]]
        a = classify(CircuitDAG.from_dict(data), self.lib, enable_overlap=True)
        data["nodes"].reverse(); data["edges"].reverse()
        b = classify(CircuitDAG.from_dict(data), self.lib, enable_overlap=True)
        self.assertEqual(a["rewrite_candidates"], b["rewrite_candidates"])
        for c in a["rewrite_candidates"]:
            self.assert_sound(c)

    def test_limits_are_partial_and_keep_sound_candidates(self):
        dag = CircuitDAG.from_dict(dag_data(2, [(op, [q]) for q in range(2) for op in ("z", "h", "x")]))
        result = classify(dag, self.lib, enable_overlap=True, overlap_limits=OverlapLimits(max_candidates=1))
        self.assertEqual(result["overlap"]["stop_reason"], "candidates_limit")
        self.assertFalse(result["overlap"]["complete_for_bounds"])
        self.assertEqual(len(result["rewrite_candidates"]), 1)
        self.assert_sound(result["rewrite_candidates"][0])
        for limit, reason in ((OverlapLimits(time_limit_seconds=1e-12), "time_limit"),
                              (OverlapLimits(max_match_attempts=1), "match_attempts_limit"),
                              (OverlapLimits(max_graph_steps=1), "graph_steps_limit"),
                              (OverlapLimits(max_views=1), "views_limit")):
            result = classify(dag, self.lib, enable_overlap=True, overlap_limits=limit)
            self.assertEqual(result["overlap"]["stop_reason"], reason)
            self.assertFalse(result["overlap"]["complete_for_bounds"])
        result = classify(self.dag, self.lib, enable_overlap=True, overlap_limits=OverlapLimits(max_original_gates=2))
        self.assertEqual(result["overlap"]["stop_reason"], "search_pruned")
        self.assertFalse(result["rewrite_candidates"])

    def test_invalid_limits_and_incomplete_base(self):
        for kwargs in ({"max_qubits": 7}, {"max_members_per_class": 0}, {"max_views": True},
                       {"time_limit_seconds": float("nan")}, {"max_virtual_gates": 257}):
            with self.subTest(kwargs=kwargs), self.assertRaises(InputError):
                classify(self.dag, self.lib, enable_overlap=True, overlap_limits=OverlapLimits(**kwargs))
        result = classify(self.dag, self.lib, Limits(max_match_attempts=1), enable_overlap=True)
        self.assertFalse(result["complete_for_bounds"])
        self.assertFalse(result["overlap"]["complete_for_bounds"])

    def test_random_context_candidates_preserve_entire_window_unitary(self):
        rng = random.Random(41)
        checked = 0
        for _ in range(16):
            steps = [("h", [0]), ("cx", [0, 1])]
            q = rng.randrange(2)
            steps.extend((op, [q]) for op in ("z", "h", "x"))
            for _ in range(3):
                op = rng.choice(["h", "x", "z", "cx"])
                qs = rng.sample(range(2), 2) if op == "cx" else [rng.randrange(2)]
                steps.append((op, qs))
            result = classify(CircuitDAG.from_dict(dag_data(2, steps)), self.lib, enable_overlap=True)
            for c in result["rewrite_candidates"]:
                self.assert_sound(c)
                checked += 1
        self.assertGreater(checked, 16)

    def test_cli_opt_in_and_exit_status(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            dag, lib = folder / "dag.json", folder / "lib.json"
            dag.write_text(json.dumps(self.dag.to_dict()))
            lib.write_text(json.dumps(self.lib.to_dict()))
            command = [sys.executable, "-m", "prefix_classifier", "--dag", str(dag), "--classes", str(lib)]
            plain = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(plain.returncode, 0, plain.stderr)
            self.assertNotIn("rewrite_candidates", json.loads(plain.stdout))
            advanced = subprocess.run(command + ["--enable-overlap"], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(advanced.returncode, 0, advanced.stderr)
            self.assertTrue(json.loads(advanced.stdout)["rewrite_candidates"])
            partial = subprocess.run(command + ["--enable-overlap", "--overlap-max-original-gates", "2"],
                                     cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(partial.returncode, 2, partial.stderr)
            self.assertEqual(json.loads(partial.stdout)["overlap"]["stop_reason"], "search_pruned")
            invalid = subprocess.run(command + ["--overlap-max-views", "2"], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(invalid.returncode, 1)
            self.assertIn("require --enable-overlap", invalid.stderr)


if __name__ == "__main__":
    unittest.main()
