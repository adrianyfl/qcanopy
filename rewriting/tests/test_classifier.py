import cmath
import copy
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from prefix_classifier import CircuitDAG, ClassLibrary, InputError, Limits, classify
from prefix_classifier.model import read_json


def dag_data(n, steps):
    """Test fixture builder only; not a QASM/front-end API."""
    nodes, edges, last = [], [], {}
    for gid, op, qs, *parameters in steps:
        nodes.append({"id": gid, "op": op, "qubits": qs,
                      "params": parameters[0] if parameters else []})
        for q in qs:
            if q in last:
                edges.append({"source": last[q], "target": gid, "qubit": q})
            last[q] = gid
    return {"schema": "qdag/1", "num_qubits": n, "nodes": nodes, "edges": edges}


def library_data():
    def member(mid, n, steps):
        return {"id": mid, "dag": dag_data(n, steps)}
    return {"schema": "qclass-library/1", "equivalence": "exact_unitary", "classes": [
        {"id": "X", "num_qubits": 1, "members": [
            member("x", 1, [("x", "x", [0])]),
            member("hzh", 1, [("h1", "h", [0]), ("z", "z", [0]), ("h2", "h", [0])])]},
        {"id": "H", "num_qubits": 1, "members": [member("h", 1, [("h", "h", [0])])]},
        {"id": "Z", "num_qubits": 1, "members": [member("z", 1, [("z", "z", [0])])]},
        {"id": "CX", "num_qubits": 2, "members": [member("cx", 2, [("cx", "cx", [0, 1])])]},
        {"id": "RZ", "num_qubits": 1, "parameters": ["theta"], "members": [
            member("rz", 1, [("rz", "rz", [0], [{"var": "theta"}])])]}]}


def run(steps, n=2, library=None, **limits):
    return classify(CircuitDAG.from_dict(dag_data(n, steps)),
                    ClassLibrary.from_dict(library or library_data()), Limits(**limits))


def atom_classes(group):
    return [a["class_id"] for a in group["signature"]["atoms"]]


class MatchingTests(unittest.TestCase):
    def test_different_members_share_fixed_class_and_composition(self):
        data = dag_data(2, [("a0", "h", [0]), ("a1", "z", [0]), ("a2", "h", [0]),
                            ("a3", "h", [0]), ("a4", "z", [0]),
                            ("b0", "x", [1]), ("b1", "h", [1]), ("b2", "z", [1])])
        dag, library = CircuitDAG.from_dict(data), ClassLibrary.from_dict(library_data())
        before = (dag.to_dict(), library.to_dict())
        rounds = []
        result = classify(dag, library, on_round=rounds.append)
        self.assertTrue(result["complete_for_bounds"])
        self.assertEqual(before, (dag.to_dict(), library.to_dict()))
        self.assertEqual(result["catalog_invariant"], {"classes_before": 5, "classes_after": 5,
                                                      "members_before": 6, "members_after": 6})
        matches = [m for m in result["matches"] if m["class_id"] == "X" and m["boundary"]["prefix"]]
        self.assertEqual(len(matches), 2)
        group = next(g for g in result["composition_groups"] if atom_classes(g) == ["X", "H", "Z"])
        self.assertEqual(len(group["region_ids"]), 2)
        regions = {r["id"]: r for r in result["regions"]}
        self.assertEqual({len(regions[r]["gate_ids"]) for r in group["region_ids"]}, {3, 5})
        self.assertTrue(any(row["round"] >= 2 for row in rounds))
        self.assertTrue(all(r["parent_region"] in regions for r in result["regions"] if r["parent_region"]))

    def test_parameters_bind_exactly(self):
        result = run([("a", "rz", [0], [0.25]), ("b", "rz", [1], [0.5]),
                      ("c", "rz", [2], [0.25])], n=3)
        groups = [g for g in result["composition_groups"] if atom_classes(g) == ["RZ"]]
        self.assertEqual(sorted(len(g["region_ids"]) for g in groups), [1, 2])
        result = run([("a", "rz", [0], [0.1]), ("b", "rz", [1], [0.1000000000000001])])
        self.assertEqual(len(result["composition_groups"]), 2)

    def test_repeated_variable_must_have_same_binding(self):
        lib = {"schema": "qclass-library/1", "equivalence": "exact_unitary", "classes": [
            {"id": "twice", "num_qubits": 1, "parameters": ["t"], "members": [
                {"id": "pair", "dag": dag_data(1, [("a", "rz", [0], [{"var": "t"}]),
                                                     ("b", "rz", [0], [{"var": "t"}])])}]}]}
        result = run([("a", "rz", [0], [1]), ("b", "rz", [0], [2]),
                      ("c", "rz", [1], [1]), ("d", "rz", [1], [1.0])], library=lib)
        self.assertEqual(len(result["matches"]), 1)
        self.assertEqual(result["matches"][0]["gate_ids"], ["c", "d"])

    def test_control_target_port_mapping(self):
        result = run([("a", "x", [0]), ("b", "cx", [0, 1]),
                      ("c", "x", [3]), ("d", "cx", [2, 3])], n=4)
        groups = [g for g in result["composition_groups"] if atom_classes(g) == ["X", "CX"]]
        self.assertEqual(len(groups), 2)
        self.assertTrue(all(len(g["region_ids"]) == 1 for g in groups))
        mappings = {tuple(m["gate_ids"]): m["qubit_map"] for m in result["matches"] if m["class_id"] == "CX"}
        self.assertEqual(mappings[("b",)], [0, 1])
        self.assertEqual(mappings[("d",)], [2, 3])

    def test_wire_renaming_and_parallel_order(self):
        # Different independent-node order and wire numbers, identical macro DAG.
        result = run([("z1", "h", [0]), ("a1", "x", [1]), ("c1", "cx", [0, 1]),
                      ("a2", "h", [3]), ("z2", "x", [2]), ("c2", "cx", [3, 2])], n=4)
        complete = [r for r in result["regions"] if len(r["gate_ids"]) == 3]
        self.assertTrue(complete)
        signatures = {r["composition_key"] for r in complete}
        self.assertEqual(len(signatures), 1)
        self.assertEqual({frozenset(r["gate_ids"]) for r in complete},
                         {frozenset(["z1", "a1", "c1"]), frozenset(["a2", "z2", "c2"])})

    def test_do_not_skip_unmatched_gate_inside_pattern(self):
        result = run([("a", "h", [0]), ("u", "unknown", [0]), ("b", "z", [0]), ("c", "h", [0])], n=1)
        self.assertFalse(any(m["class_id"] == "X" for m in result["matches"]))
        self.assertIn("u", result["unmatched_gate_ids"])
        self.assertFalse(any({"a", "b"}.issubset(r["gate_ids"]) for r in result["regions"]))

    def test_reject_nonconvex_pattern_match(self):
        lib = {"schema": "qclass-library/1", "equivalence": "exact_unitary", "classes": [
            {"id": "HX", "num_qubits": 2, "members": [
                {"id": "parallel", "dag": dag_data(2, [("h", "h", [0]), ("x", "x", [1])])}]}]}
        # H -> CX(0,2) -> CX(2,1) -> X: no missing same-wire internal edge,
        # but the candidate {H,X} has a path that leaves and reenters it.
        result = run([("h", "h", [0]), ("u", "cx", [0, 2]),
                      ("v", "cx", [2, 1]), ("x", "x", [1])], n=3, library=lib)
        self.assertEqual(result["matches"], [])

    def test_cross_wire_predecessor_prevents_false_prefix(self):
        result = run([("u", "unknown", [1]), ("a", "h", [0]), ("b", "cx", [0, 1]),
                      ("v", "unknown", [0])])
        candidates = [r for r in result["regions"] if set(r["gate_ids"]) == {"a", "b"}]
        self.assertTrue(candidates)
        self.assertTrue(all(not r["boundary"]["prefix"] for r in candidates))
        cx = next(m for m in result["matches"] if m["class_id"] == "CX")
        self.assertFalse(cx["boundary"]["prefix"])

    def test_suffix_growth_keeps_forward_operator_order(self):
        result = run([("u", "unknown", [0]), ("a", "x", [0]),
                      ("b", "h", [0]), ("c", "z", [0])], n=1)
        region = next(r for r in result["regions"] if set(r["gate_ids"]) == {"a", "b", "c"})
        self.assertEqual(region["boundary"], {"prefix": False, "suffix": True})
        group = next(g for g in result["composition_groups"] if g["key"] == region["composition_key"])
        self.assertEqual(atom_classes(group), ["X", "H", "Z"])

    def test_fan_in_and_fan_out_do_not_require_premature_boundary_closure(self):
        result = run([("h0", "h", [0]), ("x1", "x", [1]), ("cx", "cx", [0, 1]),
                      ("z0", "z", [0]), ("h1", "h", [1])])
        by_gates = {frozenset(r["gate_ids"]): r for r in result["regions"]}
        self.assertIn(frozenset(["h0", "x1", "cx"]), by_gates)
        self.assertTrue(by_gates[frozenset(["h0", "x1", "cx"])]["boundary"]["prefix"])
        self.assertIn(frozenset(["cx", "z0", "h1"]), by_gates)
        self.assertTrue(by_gates[frozenset(["cx", "z0", "h1"])]["boundary"]["suffix"])
        self.assertIn(frozenset(["h0", "x1", "cx", "z0", "h1"]), by_gates)
        partial = by_gates[frozenset(["h0", "cx"])]
        self.assertEqual(partial["boundary"], {"prefix": False, "suffix": False})

    def test_link_witnesses_cannot_be_stitched_across_positions(self):
        result = run([("a", "x", [0]), ("b", "h", [0]),
                      ("c", "h", [1]), ("d", "z", [1])])
        edges = {(l["source_class"], l["target_class"]) for l in result["links"]}
        self.assertIn(("X", "H"), edges)
        self.assertIn(("H", "Z"), edges)
        self.assertFalse(any(atom_classes(g) == ["X", "H", "Z"] for g in result["composition_groups"]))
        self.assertTrue(all(len(r["qubit_map"]) == 1 for r in result["regions"]))

    def test_repeated_class_edges_terminate_and_keep_children(self):
        result = run([(str(i), "h", [0]) for i in range(6)], n=1)
        self.assertTrue(result["complete_for_bounds"])
        self.assertEqual(result["catalog_invariant"]["classes_after"], 5)
        self.assertTrue(any(l["source_class"] == l["target_class"] == "H" for l in result["links"]))
        self.assertEqual(max(len(r["gate_ids"]) for r in result["regions"]), 6)
        self.assertTrue(any(len(r["gate_ids"]) == 1 for r in result["regions"]))
        for r in result["regions"]:
            self.assertEqual(sum(len(next(m for m in result["matches"] if m["id"] == mid)["gate_ids"])
                                 for mid in r["parts"]), len(r["gate_ids"]))

    def test_input_array_order_does_not_change_classification(self):
        data = dag_data(2, [("a", "x", [0]), ("b", "h", [0]), ("c", "cx", [0, 1])])
        lib = ClassLibrary.from_dict(library_data())
        first = classify(CircuitDAG.from_dict(data), lib)
        data["nodes"].reverse()
        data["edges"].reverse()
        second = classify(CircuitDAG.from_dict(data), lib)
        for key in ("matches", "regions", "composition_groups", "links", "dag_fingerprint"):
            self.assertEqual(first[key], second[key])

    def test_empty_circuit_or_catalog(self):
        dag = CircuitDAG.from_dict(dag_data(2, []))
        lib = ClassLibrary.from_dict(library_data())
        result = classify(dag, lib)
        self.assertTrue(result["complete_for_bounds"])
        self.assertEqual(result["matches"], [])
        lib = ClassLibrary.from_dict({"schema": "qclass-library/1", "equivalence": "exact_unitary", "classes": []})
        result = classify(CircuitDAG.from_dict(dag_data(1, [("a", "h", [0])])), lib)
        self.assertEqual(result["unmatched_gate_ids"], ["a"])
        self.assertEqual(result["catalog_invariant"]["classes_after"], 0)

    def test_global_phase_is_preserved(self):
        data = dag_data(1, [("a", "h", [0])])
        data["global_phase"] = "pi/7"
        self.assertEqual(classify(CircuitDAG.from_dict(data), ClassLibrary.from_dict(library_data()))["global_phase"], "pi/7")

    def test_composition_groups_have_equal_unitary_actions(self):
        steps = [("a0", "h", [0]), ("a1", "z", [0]), ("a2", "h", [0]),
                 ("a3", "h", [0]), ("a4", "cx", [0, 1]), ("a5", "rz", [1], [0.3]),
                 ("b0", "x", [3]), ("b1", "h", [3]), ("b2", "cx", [3, 2]), ("b3", "rz", [2], [0.3])]
        dag = CircuitDAG.from_dict(dag_data(4, steps))
        result = classify(dag, ClassLibrary.from_dict(library_data()))
        regions = {r["id"]: r for r in result["regions"]}

        def simulate(region, basis):
            wires = region["qubit_map"]
            mapping = {q: i for i, q in enumerate(wires)}
            state = [0j] * (1 << len(wires))
            state[basis] = 1 + 0j
            for gid in region["gate_ids"]:
                gate = dag.by_id[gid]
                q = mapping[gate.qubits[0]]
                if gate.op == "cx":
                    target = mapping[gate.qubits[1]]
                    out = [0j] * len(state)
                    for i, amplitude in enumerate(state):
                        out[i ^ (1 << target) if i & (1 << q) else i] = amplitude
                    state = out
                    continue
                if gate.op == "h":
                    h = 1 / math.sqrt(2)
                    matrix = [[h, h], [h, -h]]
                elif gate.op == "x":
                    matrix = [[0, 1], [1, 0]]
                elif gate.op == "z":
                    matrix = [[1, 0], [0, -1]]
                else:
                    t = gate.params[0]
                    matrix = [[cmath.exp(-0.5j * t), 0], [0, cmath.exp(0.5j * t)]]
                out = [0j] * len(state)
                for i in range(len(state)):
                    row = (i >> q) & 1
                    for col in range(2):
                        source = (i & ~(1 << q)) | (col << q)
                        out[i] += matrix[row][col] * state[source]
                state = out
            return state

        compared = 0
        for group in result["composition_groups"]:
            if len(group["region_ids"]) < 2:
                continue
            for basis in range(1 << group["signature"]["num_qubits"]):
                states = [simulate(regions[r], basis) for r in group["region_ids"]]
                for state in states[1:]:
                    self.assertLess(max(abs(a - b) for a, b in zip(states[0], state)), 1e-12)
                    compared += 1
        self.assertGreater(compared, 5)


class ValidationAndBudgetTests(unittest.TestCase):
    def test_malformed_dags(self):
        base = dag_data(2, [("a", "h", [0]), ("b", "cx", [0, 1]), ("c", "z", [1])])
        bad = []
        missing = copy.deepcopy(base); missing["edges"].pop(); bad.append(missing)
        duplicate = copy.deepcopy(base); duplicate["edges"].append(duplicate["edges"][0]); bad.append(duplicate)
        endpoint = copy.deepcopy(base); endpoint["edges"][0]["source"] = "missing"; bad.append(endpoint)
        operand = copy.deepcopy(base); operand["nodes"][1]["qubits"] = [0, 0]; bad.append(operand)
        unknown = copy.deepcopy(base); unknown["nodes"][0]["condition"] = "c0"; bad.append(unknown)
        nonunitary = copy.deepcopy(base); nonunitary["nodes"][0]["op"] = "measure"; bad.append(nonunitary)
        cycle = dag_data(2, [("a", "cx", [0, 1]), ("b", "cx", [0, 1])])
        cycle["edges"][1] = {"source": "b", "target": "a", "qubit": 1}; bad.append(cycle)
        for data in bad:
            with self.subTest(data=data), self.assertRaises(InputError):
                CircuitDAG.from_dict(data)

    def test_invalid_catalogs(self):
        base = library_data()
        bad = []
        phase = copy.deepcopy(base); phase["equivalence"] = "approximate"; bad.append(phase)
        unknown = copy.deepcopy(base); unknown["classes"][-1]["parameters"] = []; bad.append(unknown)
        duplicate = copy.deepcopy(base); duplicate["classes"].append(duplicate["classes"][0]); bad.append(duplicate)
        empty = copy.deepcopy(base); empty["classes"][0]["members"] = []; bad.append(empty)
        for data in bad:
            with self.subTest(data=data), self.assertRaises(InputError):
                ClassLibrary.from_dict(data)

    def test_limits_return_sound_partial_results(self):
        steps = [(str(i), "h", [0]) for i in range(8)]
        cases = [(dict(max_rounds=0), "rounds_limit"),
                 (dict(max_matches=1), "matches_limit"),
                 (dict(max_regions=1), "regions_limit"),
                 (dict(max_match_attempts=1), "match_attempts_limit"),
                 (dict(max_extension_attempts=1), "extension_attempts_limit"),
                 (dict(max_link_witnesses=1), "link_witnesses_limit"),
                 (dict(time_limit_seconds=1e-12), "time_limit")]
        for limits, expected in cases:
            with self.subTest(limits=limits):
                result = run(steps, n=1, **limits)
                self.assertEqual(result["stop_reason"], expected)
                self.assertFalse(result["complete_for_bounds"])
                match_ids = {m["id"] for m in result["matches"]}
                region_ids = {r["id"] for r in result["regions"]}
                for region in result["regions"]:
                    self.assertTrue(set(region["parts"]) <= match_ids)
                    self.assertTrue(region["parent_region"] is None or region["parent_region"] in region_ids)
                for link in result["links"]:
                    for witness in link["witnesses"]:
                        self.assertIn(witness["source_match"], match_ids)
                        self.assertIn(witness["target_match"], match_ids)

    def test_fixed_point_is_relative_to_size_bounds(self):
        result = run([(str(i), "h", [0]) for i in range(5)], n=1, max_gates=2)
        self.assertTrue(result["complete_for_bounds"])
        self.assertGreater(result["statistics"]["size_pruned"], 0)
        self.assertLessEqual(max(len(r["gate_ids"]) for r in result["regions"]), 2)

    def test_duplicate_json_keys_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "bad.json"
            p.write_text('{"schema":"qdag/1","schema":"other"}')
            with self.assertRaises(InputError):
                read_json(p)

    def test_cli_end_to_end_and_input_overwrite_guard(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            dag, lib, output = folder / "dag.json", folder / "classes.json", folder / "result.json"
            dag.write_text(json.dumps(dag_data(1, [("a", "x", [0]), ("b", "h", [0])])))
            lib.write_text(json.dumps(library_data()))
            command = [sys.executable, "-m", "prefix_classifier", "--dag", str(dag), "--classes", str(lib)]
            done = subprocess.run(command + ["-o", str(output), "--progress"], cwd=root, capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(output.read_text())["schema"], "qclassification/1")
            self.assertTrue(done.stderr)
            partial = subprocess.run(command + ["--max-rounds", "0"], cwd=root, capture_output=True, text=True)
            self.assertEqual(partial.returncode, 2)
            self.assertFalse(json.loads(partial.stdout)["complete_for_bounds"])
            before = dag.read_text()
            failed = subprocess.run(command + ["-o", str(dag)], cwd=root, capture_output=True, text=True)
            self.assertEqual(failed.returncode, 1)
            self.assertEqual(dag.read_text(), before)


if __name__ == "__main__":
    unittest.main()
