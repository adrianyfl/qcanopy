import cmath
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from prefix_classifier import (CircuitDAG, InputError, classify, load_dag,
                               load_library, remove_swaps)


ROOT = Path(__file__).resolve().parents[1]


def make_dag(n, steps, phase=0):
    nodes, edges, last = [], [], {}
    for i, (op, qs, params) in enumerate(steps):
        gid = f"g{i}"
        nodes.append({"id": gid, "op": op, "qubits": qs, "params": params})
        for q in qs:
            if q in last:
                edges.append({"source": last[q], "target": gid, "qubit": q})
            last[q] = gid
    return CircuitDAG.from_dict({"schema": "qdag/1", "num_qubits": n,
                                 "nodes": nodes, "edges": edges, "global_phase": phase})


def simulate(steps, initial, phase):
    """Independent tiny reference: SWAP physically permutes amplitudes."""
    state = list(initial)
    for op, qs, params in steps:
        out = [0j] * len(state)
        if op in ("swap", "cx"):
            a, b = qs
            for i, value in enumerate(state):
                if op == "swap":
                    j = i ^ ((1 << a) | (1 << b)) if ((i >> a) ^ (i >> b)) & 1 else i
                else:
                    j = i ^ (1 << b) if (i >> a) & 1 else i
                out[j] += value
        else:
            q = qs[0]
            if op == "h":
                t = 1 / math.sqrt(2)
                matrix = ((t, t), (t, -t))
            elif op == "x":
                matrix = ((0, 1), (1, 0))
            elif op == "rz":
                t = params[0]
                matrix = ((cmath.exp(-0.5j * t), 0), (0, cmath.exp(0.5j * t)))
            else:
                raise AssertionError(f"unsupported reference gate: {op}")
            for i, value in enumerate(state):
                col = (i >> q) & 1
                for row in range(2):
                    j = (i & ~(1 << q)) | (row << q)
                    out[j] += matrix[row][col] * value
        state = out
    return [cmath.exp(1j * phase) * value for value in state]


class SwapRemovalTests(unittest.TestCase):
    def test_rebuild_edges_and_preserve_input(self):
        dag = make_dag(2, [("h", [0], []), ("swap", [0, 1], []),
                           ("z", [0], []), ("x", [1], [])])
        before = dag.to_dict()
        result = remove_swaps(dag)
        self.assertEqual(dag.to_dict(), before)
        self.assertEqual(result.output_wire_to_core, (1, 0))
        self.assertEqual(result.removed_gate_ids, ("g1",))
        self.assertEqual(result.core_dag.by_id["g2"].qubits, (1,))
        self.assertEqual(result.core_dag.by_id["g3"].qubits, (0,))
        self.assertEqual(result.core_dag.edges, (("g0", "g3", 0),))
        self.assertEqual(set(result.core_dag.by_id), {"g0", "g2", "g3"})

    def test_three_cycle_mapping_direction_and_ordered_operands(self):
        result = remove_swaps(make_dag(3, [("swap", [0, 1], []),
                                           ("swap", [1, 2], []), ("cx", [0, 2], [])]))
        self.assertEqual(result.output_wire_to_core, (1, 2, 0))
        self.assertEqual(result.core_dag.by_id["g2"].qubits, (1, 0))
        self.assertEqual([result.output_index_to_core(i) for i in range(8)],
                         [0, 2, 4, 6, 1, 3, 5, 7])
        # Core bits are indexed by qubit, not by printed big-endian bitstring order.
        core_bits = [1, 0, 1]
        self.assertEqual([core_bits[q] for q in result.output_wire_to_core], [0, 1, 1])

    def test_no_swaps_empty_and_swap_only(self):
        for steps in ([], [("h", [0], []), ("cx", [0, 1], [])]):
            dag = make_dag(3, steps, phase=0.25)
            result = remove_swaps(dag)
            self.assertEqual(result.core_dag.to_dict(), dag.to_dict())
            self.assertEqual(result.output_wire_to_core, (0, 1, 2))
            self.assertFalse(result.removed_gate_ids)
        result = remove_swaps(make_dag(4, [("swap", [2, 0], []), ("swap", [0, 1], [])]))
        self.assertEqual(result.core_dag.num_qubits, 4)
        self.assertFalse(result.core_dag.nodes)
        self.assertEqual(result.output_wire_to_core, (1, 2, 0, 3))

    def test_preserve_phase_parameters_and_other_gate_semantics(self):
        dag = make_dag(3, [("swap", [0, 1], []),
                           ("custom", [0, 2], ["theta", 0.125]),
                           ("iswap", [0, 1], []), ("cswap", [2, 0, 1], []),
                           ("SWAP", [0, 1], [])], phase="phi")
        result = remove_swaps(dag)
        self.assertEqual(result.core_dag.global_phase, "phi")
        for gid in ("g1", "g2", "g3", "g4"):
            self.assertEqual(result.core_dag.by_id[gid].op, dag.by_id[gid].op)
            self.assertEqual(result.core_dag.by_id[gid].params, dag.by_id[gid].params)
        self.assertEqual(result.core_dag.by_id["g1"].qubits, (1, 2))
        self.assertEqual(result.core_dag.by_id["g3"].qubits, (2, 1, 0))
        triple = make_dag(2, [("cx", [0, 1], []), ("cx", [1, 0], []), ("cx", [0, 1], [])])
        self.assertEqual(remove_swaps(triple).core_dag.to_dict(), triple.to_dict())

    def test_input_array_order_is_irrelevant(self):
        dag = make_dag(4, [("h", [0], []), ("x", [3], []), ("swap", [0, 2], []),
                           ("swap", [1, 3], []), ("cx", [3, 2], [])])
        data = dag.to_dict()
        data["nodes"].reverse()
        data["edges"].reverse()
        self.assertEqual(remove_swaps(dag).to_dict(),
                         remove_swaps(CircuitDAG.from_dict(data)).to_dict())

    def test_reject_malformed_explicit_swap(self):
        for qs, params in (([0], []), ([0, 1, 2], []), ([0, 1], [0])):
            with self.subTest(qs=qs, params=params), self.assertRaises(InputError):
                remove_swaps(make_dag(3, [("swap", qs, params)]))

    def test_index_bounds(self):
        result = remove_swaps(make_dag(2, []))
        for index in (-1, 4, 1.0, True):
            with self.subTest(index=index), self.assertRaises(InputError):
                result.output_index_to_core(index)

    def test_exact_unitary_action_on_all_basis_and_complex_states(self):
        rng = random.Random(597)
        for n in (2, 3, 4):
            for case in range(12):
                # Force entanglement before a SWAP, then exercise arbitrary interactions.
                steps = [("h", [0], []), ("cx", [0, 1], []), ("swap", [0, 1], [])]
                for _ in range(17):
                    op = rng.choice(["h", "x", "rz", "cx", "swap"])
                    qs = rng.sample(range(n), 2) if op in ("cx", "swap") else [rng.randrange(n)]
                    steps.append((op, qs, [rng.uniform(-3, 3)] if op == "rz" else []))
                phase = rng.uniform(-3, 3)
                result = remove_swaps(make_dag(n, steps, phase))
                core = result.core_dag
                core_steps = [(core.by_id[g].op, core.by_id[g].qubits, core.by_id[g].params)
                              for g in core.order]
                states = [[complex(i == basis) for i in range(1 << n)] for basis in range(1 << n)]
                values = [complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(1 << n)]
                norm = math.sqrt(sum(abs(v) ** 2 for v in values))
                states.append([v / norm for v in values])
                for initial in states:
                    expected = simulate(steps, initial, phase)
                    actual = simulate(core_steps, initial, core.global_phase)
                    restored = [actual[result.output_index_to_core(i)] for i in range(1 << n)]
                    with self.subTest(n=n, case=case):
                        self.assertLess(max(abs(a - b) for a, b in zip(expected, restored)), 1e-12)

    def test_serialization_and_downstream_classification(self):
        prepared = remove_swaps(load_dag(ROOT / "examples/swaps.json"))
        encoded = json.loads(json.dumps(prepared.to_dict()))
        self.assertEqual(encoded["schema"], "qswap-removal/1")
        self.assertEqual(encoded["output_wire_to_core"], [1, 2, 0])
        core = CircuitDAG.from_dict(encoded["core_dag"])
        result = classify(core, load_library(ROOT / "examples/classes.json"))
        self.assertTrue(result["complete_for_bounds"])
        self.assertEqual(result["unmatched_gate_ids"], [])
        tail = next(m for m in result["matches"] if m["gate_ids"] == ["tail"])
        self.assertEqual(tail["qubit_map"], [1, 0])
        second = remove_swaps(core)
        self.assertEqual(second.core_dag.to_dict(), core.to_dict())
        self.assertEqual(second.output_wire_to_core, (0, 1, 2))

    def test_cli_roundtrip_and_error_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            source, output = folder / "input.json", folder / "prepared.json"
            source.write_text((ROOT / "examples/swaps.json").read_text())
            command = [sys.executable, "-m", "prefix_classifier.swap_cli", "--dag", str(source)]
            completed = subprocess.run(command + ["-o", str(output)], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            expected = remove_swaps(load_dag(source)).to_dict()
            self.assertEqual(json.loads(output.read_text()), expected)
            stdout = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(stdout.returncode, 0, stdout.stderr)
            self.assertEqual(json.loads(stdout.stdout), expected)
            before = source.read_text()
            failed = subprocess.run(command + ["-o", str(source)], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(failed.returncode, 1)
            self.assertEqual(source.read_text(), before)
            source.write_text(json.dumps(make_dag(2, [("swap", [0], [])]).to_dict()))
            output.unlink()
            failed = subprocess.run(command + ["-o", str(output)], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(failed.returncode, 1)
            self.assertIn("swap requires", failed.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
