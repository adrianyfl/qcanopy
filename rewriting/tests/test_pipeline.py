"""Integration checks for the public team-facing entry points."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from prefix_classifier import ClassLibrary, CircuitDAG, load_dag, load_library
from prefix_classifier.pipeline import analyze_dag, compare_dags
from prefix_classifier.qcanopy import sequence_to_dag


ROOT = Path(__file__).resolve().parents[1]


def dag(n, *gates):
    return sequence_to_dag(n, [dict(name=op, qubits=list(qs), params=[]) for op, qs in gates])


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.library = load_library(ROOT / "examples/classes.json")

    def test_single_analysis_preserves_mapping_and_classification_coordinates(self):
        source = load_dag(ROOT / "examples/swaps.json")
        before = source.to_dict()
        result = analyze_dag(source, self.library)
        prepared, c = result["swap_removal"], result["classification"]
        self.assertEqual(prepared["output_wire_to_core"], [1, 2, 0])
        core_ids = {g["id"] for g in prepared["core_dag"]["nodes"]}
        self.assertEqual(set(c["annotations"]), core_ids)
        self.assertFalse(any(g["op"] == "swap" for g in prepared["core_dag"]["nodes"]))
        self.assertNotIn("rewrite_candidates", c)
        self.assertEqual(before, source.to_dict())

    def test_overlap_is_opt_in_and_proposals_remain_separate(self):
        source = load_dag(ROOT / "examples/overlap-circuit.json")
        lib = load_library(ROOT / "examples/overlap-classes.json")
        off = analyze_dag(source, lib)
        on = analyze_dag(source, lib, enable_overlap=True)
        self.assertNotIn("rewrite_candidates", off["classification"])
        self.assertTrue(on["classification"]["rewrite_candidates"])
        self.assertEqual(off["swap_removal"], on["swap_removal"])

    def test_pair_layout_mismatch_cannot_erase_a_real_swap(self):
        lib = ClassLibrary.from_dict(dict(schema="qclass-library/1", equivalence="exact_unitary", classes=[]))
        source = dag(2, ("swap", (0, 1)))
        r = compare_dags(dag(2), source, lib)
        self.assertFalse(r["layout"]["cores_used"])
        self.assertEqual(r["proof"]["verdict"], "inconclusive")
        self.assertTrue(r["proof"]["residual_right"]["nodes"])
        r = compare_dags(source, source, lib)
        self.assertTrue(r["layout"]["cores_used"])
        self.assertEqual(r["proof"]["verdict"], "equivalent_up_to_global_phase")

    def test_pair_inverse_rule_without_overlap(self):
        lib = ClassLibrary.from_dict(dict(schema="qclass-library/1", equivalence="exact_unitary", classes=[]))
        r = compare_dags(dag(1), dag(1, ("s", (0,)), ("sdg", (0,))), lib)
        self.assertEqual(r["proof"]["verdict"], "equivalent_up_to_global_phase")
        self.assertTrue(any(p["kind"] == "inverse_identity" for p in r["proof"]["proof"]))

    def test_standalone_cli_ir_bridge_and_error_codes(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            ir = folder / "source.ir.json"
            ir.write_text(json.dumps(dict(schema_version=1, n=2, global_phase=.25, gates=[
                dict(name="swap", qubits=[0, 1], params=[])])))
            cmd = [sys.executable, str(ROOT / "run.py")]
            r = subprocess.run(cmd + ["remove-swaps", "--dag", str(ir), "--input-format", "qcanopy-ir"],
                               cwd=folder, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            prepared = json.loads(r.stdout)
            self.assertEqual(prepared["output_wire_to_core"], [1, 0])
            self.assertEqual(prepared["core_dag"]["global_phase"], .25)
            empty = folder / "empty.json"
            empty.write_text(json.dumps(dag(2).to_dict()))
            source = folder / "swap.json"
            source.write_text(json.dumps(dag(2, ("swap", (0, 1))).to_dict()))
            out = folder / "results/pair.json"
            r = subprocess.run(cmd + ["pair", "--left", str(empty), "--right", str(source), "-o", str(out)],
                               cwd=folder, capture_output=True, text=True)
            self.assertEqual(r.returncode, 2, r.stderr)
            self.assertEqual(json.loads(out.read_text())["proof"]["verdict"], "inconclusive")
            before = source.read_bytes()
            r = subprocess.run(cmd + ["classify", "--dag", str(source), "-o", str(source)],
                               cwd=folder, capture_output=True, text=True)
            self.assertEqual(r.returncode, 1)
            self.assertEqual(source.read_bytes(), before)
            r = subprocess.run(cmd + ["classify", "--dag", str(source), "--time-limit-seconds", "nan"],
                               cwd=folder, capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)


if __name__ == "__main__":
    unittest.main()
