import unittest
from pathlib import Path

from prefix_classifier import load_library, Limits, OverlapLimits
from prefix_classifier.pair import prove_pair, structural_peel
from prefix_classifier.qcanopy import sequence_to_dag, parse_qcanopy_ir
from prefix_classifier.model import InputError


def d(n,*ops):
    return sequence_to_dag(n,[dict(name=o,qubits=list(q),params=[]) for o,q in ops])


class PairProofTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib=load_library(Path(__file__).resolve().parents[1]/'rules/verified/classes.phase.json')

    def test_independent_order_is_already_structurally_identical(self):
        a=d(2,('h',[0]),('x',[1]),('cx',[0,1]))
        b=d(2,('x',[1]),('h',[0]),('cx',[0,1]))
        p=structural_peel(a,b)
        self.assertTrue(p['equivalent']);self.assertEqual(p['removed'],6)

    def test_fixed_class_proof_and_identity_region(self):
        r=prove_pair(d(1,('s',[0])),d(1,('t',[0]),('t',[0])),self.lib)
        self.assertEqual(r['verdict'],'equivalent_up_to_global_phase')
        self.assertTrue(any(p['kind']=='library_pair' for p in r['proof']))
        r=prove_pair(d(1),d(1,('h',[0]),('h',[0])),self.lib,enable_inverse=False)
        self.assertEqual(r['verdict'],'equivalent_up_to_global_phase')
        self.assertTrue(any(p['kind']=='library_identity' for p in r['proof']))

    def test_ordered_physical_ports_cannot_be_renamed_between_circuits(self):
        r=prove_pair(d(2,('cx',[0,1])),d(2,('cx',[1,0])),self.lib)
        self.assertEqual(r['verdict'],'inconclusive')

    def test_partial_overlap_proof_without_truth_or_numeric_oracle(self):
        a=d(1,('sdg',[0]),('sx',[0]),('s',[0]),('x',[0]))
        b=d(1,('z',[0]),('h',[0]),('x',[0]))
        saved=(a.to_dict(),b.to_dict())
        lim=Limits(max_qubits=1,max_gates=16,time_limit_seconds=2)
        r=prove_pair(a,b,self.lib,lim,OverlapLimits(max_qubits=1,max_candidates=500,
            max_views=3000,max_members_per_class=8,max_views_per_match=24,max_match_attempts=1000000),
            enable_overlap=True,overlap_seconds=5)
        self.assertFalse(r['base_proved'])
        self.assertEqual(r['verdict'],'equivalent_up_to_global_phase')
        self.assertTrue(any(p['kind']=='overlap_substitution' for p in r['proof']))
        self.assertEqual(saved,(a.to_dict(),b.to_dict()))

    def test_non_equivalent_and_same_zero_state_are_not_proved(self):
        for a,b in [(d(1),d(1,('z',[0]))),(d(2),d(2,('swap',[0,1]))),
                    (d(1,('h',[0])),d(1,('x',[0])))]:
            self.assertEqual(prove_pair(a,b,self.lib)['verdict'],'inconclusive')

    def test_explicit_global_phase_is_allowed_and_residual_budget_is_not_success(self):
        a=sequence_to_dag(1,[],phase=.123);b=sequence_to_dag(1,[],phase=-.456)
        self.assertEqual(prove_pair(a,b,self.lib)['verdict'],'equivalent_up_to_global_phase')
        r=prove_pair(d(1,('h',[0])),d(1,('x',[0])),self.lib,max_residual_gates=1)
        self.assertEqual(r['verdict'],'inconclusive')

    def test_opaque_matrix_import_is_opt_in_and_keeps_exact_entries(self):
        raw=dict(schema_version=1,n=1,global_phase=0,gates=[dict(name='unitary_1q',
            qubits=[0],params=[],matrix=[[1,0],[0,0],[0,0],[0,1]])])
        with self.assertRaises(InputError):parse_qcanopy_ir(raw)
        a=parse_qcanopy_ir(raw,allow_opaque_matrices=True)
        self.assertEqual(a.nodes[0].params,(1,0,0,0,0,0,0,1))
        self.assertTrue(structural_peel(a,a)['equivalent'])
        raw['gates'][0]['matrix'].pop()
        with self.assertRaises(InputError):parse_qcanopy_ir(raw,allow_opaque_matrices=True)
