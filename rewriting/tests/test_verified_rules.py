import json
from pathlib import Path
import sys
import unittest

from prefix_classifier import ClassLibrary, CircuitDAG, classify, load_library

ROOT = Path(__file__).resolve().parents[1]
RULES = ROOT / "rules/verified"
sys.path.insert(0, str(RULES))
from exact import unitary, compare, gate_matrix
from build import dag, minimal_pair, rust_param
from verify import verify_gate


def g(op, *qs, params=()):
    return op, tuple(qs), tuple(params)


class VerifiedRuleTests(unittest.TestCase):
    def test_exact_and_projective_equalities_are_distinct(self):
        self.assertEqual(compare(unitary(1, (g("h", 0), g("h", 0))), unitary(1, (g("id", 0),))), "exact")
        self.assertEqual(compare(unitary(1, (g("t", 0), g("x", 0))),
                                 unitary(1, (g("x", 0), g("tdg", 0)))), "projective")
        self.assertEqual(compare(unitary(1, (g("rz", 0, params=("$theta",)),)),
                                 unitary(1, (g("p", 0, params=("$theta",)),))), "projective")
        self.assertEqual(compare(unitary(1, (g("id", 0),)), unitary(1, (g("z", 0),))), "different")
        self.assertEqual(compare(unitary(2, (g("cx", 0, 1),)), unitary(2, (g("cx", 1, 0),))), "different")

    def test_near_rotations_are_not_accepted_by_a_tolerance(self):
        with self.assertRaises(ValueError): gate_matrix("rx", (1e-12,))
        self.assertEqual(compare(unitary(1, (g("rx", 0, params=(0,)),)),
                                 unitary(1, (g("rx", 0, params=("pi/4",)),))), "different")
        with self.assertRaises(ValueError): rust_param("Param::ParameterExpression(theta_div_2.clone())")

    def test_symbolic_rotation_conjugation_and_unitarity(self):
        a = (g("h", 0), g("rz", 0, params=("$theta",)), g("h", 0))
        b = (g("rx", 0, params=("$theta",)),)
        self.assertEqual(compare(unitary(1, a), unitary(1, b)), "exact")
        for op, params in [("u", ("$theta", "$phi", "$lam")), ("ccx", ()),
                           ("cswap", ()), ("ch", ()), ("crz", ("$theta",)),
                           ("rxx", ("$theta",)), ("ryy", ("$theta",)), ("rzz", ("$theta",))]:
            verify_gate(op, params)

    def test_dedup_preserves_interfaces_and_removes_context(self):
        a = (g("x", 0), g("cx", 0, 1))
        b = (g("cx", 0, 1), g("x", 0), g("x", 1))
        base = minimal_pair(a, b)
        shifted_a = (g("h", 2), g("x", 2), g("cx", 2, 1), g("t", 1))
        shifted_b = (g("h", 2), g("cx", 2, 1), g("x", 1), g("x", 2), g("t", 1))
        normalized = minimal_pair(shifted_b, shifted_a)
        self.assertEqual(base[:3], normalized[:3])
        self.assertEqual(normalized[3], 2)
        self.assertEqual(compare(unitary(base[0], base[1]), unitary(base[0], base[2])), "exact")

    def test_artifacts_load_and_symbolic_members_are_nonempty(self):
        phase = load_library(RULES / "classes.phase.json")
        exact = load_library(RULES / "classes.exact.json")
        self.assertEqual(phase.equivalence, "up_to_global_phase")
        self.assertEqual(exact.equivalence, "exact_unitary")
        self.assertGreater(sum(len(c.members) for c in phase.classes), 1000)
        self.assertTrue(any(c.parameters for c in phase.classes))
        for c in phase.classes:
            self.assertEqual(len({json.dumps(m.dag.to_dict(), sort_keys=True) for m in c.members}), len(c.members))
            self.assertTrue(all(m.dag.nodes for m in c.members))

    def test_phase_mode_propagates_to_contextual_overlap(self):
        data = {"schema": "qclass-library/1", "equivalence": "up_to_global_phase", "classes": [
            {"id": "move", "num_qubits": 1, "members": [
                {"id": "tx", "dag": dag(1, (g("t", 0), g("x", 0)))},
                {"id": "xtdg", "dag": dag(1, (g("x", 0), g("tdg", 0)))}]},
            {"id": "cancel", "num_qubits": 1, "members": [
                {"id": "pair", "dag": dag(1, (g("tdg", 0), g("t", 0)))},
                {"id": "id", "dag": dag(1, (g("id", 0),))}]}]}
        lib = ClassLibrary.from_dict(data)
        circuit = CircuitDAG.from_dict(dag(1, (g("t", 0), g("x", 0), g("t", 0))))
        result = classify(circuit, lib, enable_overlap=True)
        self.assertEqual(result["equivalence"], "up_to_global_phase")
        self.assertEqual(result["overlap"]["equivalence"], "up_to_global_phase")
        self.assertTrue(result["rewrite_candidates"])
        from exact import from_dag
        self.assertTrue(any(compare(unitary(1, from_dag(c["original_dag"])),
                                    unitary(1, from_dag(c["alternative_dag"]))) == "projective"
                            for c in result["rewrite_candidates"]))
        self.assertTrue(all(c["equivalence"] == "up_to_global_phase" for c in result["rewrite_candidates"]))


if __name__ == "__main__":
    unittest.main()
