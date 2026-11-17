import itertools
import unittest

from prefix_classifier import Limits, classify, load_library
from prefix_classifier.boundary import boundary_windows
from prefix_classifier.qcanopy import sequence_to_dag
from pathlib import Path


class BoundaryWindowsTests(unittest.TestCase):
    def test_closure_and_all_small_true_boundaries_are_retained(self):
        seq=[('h',[0]),('x',[1]),('cx',[0,1]),('z',[0]),('cx',[1,2]),('h',[2]),('x',[0])]
        d=sequence_to_dag(3,[dict(name=o,qubits=q,params=[]) for o,q in seq])
        for depth in (1,2,3,4):
            ws=boundary_windows(d,depth)
            for side,w in ws.items():
                adjacency=d.predecessors if side=='prefix' else d.successors
                for g in w.order:
                    self.assertTrue(set(adjacency[g])<=set(w.order))
                self.assertEqual(set(w.edges),{e for e in d.edges if e[0] in w.by_id and e[1] in w.by_id})
            for size in range(1,depth+1):
                for subset in itertools.combinations(d.order,size):
                    pre,suf=d.boundary(subset)
                    if pre:self.assertTrue(set(subset)<=set(ws['prefix'].order))
                    if suf:self.assertTrue(set(subset)<=set(ws['suffix'].order))

    def test_extracted_tail_cut_is_not_a_true_full_circuit_suffix(self):
        d=sequence_to_dag(1,[dict(name='h',qubits=[0],params=[]) for _ in range(10)])
        ws=boundary_windows(d,2)
        self.assertEqual(set(ws['prefix'].order),{'g0','g1'})
        self.assertTrue(ws['prefix'].boundary(['g1'])[1])
        self.assertFalse(d.boundary(['g1'])[1])
        self.assertEqual(set(ws['suffix'].order),{'g8','g9'})

    def test_full_vs_window_boundary_region_keys(self):
        path=Path(__file__).resolve().parents[1]/'rules/verified/classes.phase.json'
        lib=load_library(path)
        d=sequence_to_dag(2,[dict(name=o,qubits=q,params=[]) for o,q in
            [('z',[0]),('h',[0]),('x',[0]),('h',[1]),('cx',[0,1]),('x',[1]),('h',[1])]])
        lim=Limits(max_qubits=1,max_gates=3,time_limit_seconds=10)
        full=classify(d,lib,lim)
        self.assertTrue(full['complete_for_bounds'])
        for side,w in boundary_windows(d,3).items():
            result=classify(w,lib,lim)
            expected={(frozenset(r['gate_ids']),r['composition_key']) for r in full['regions'] if r['boundary'][side]}
            actual={(frozenset(r['gate_ids']),r['composition_key']) for r in result['regions']
                    if d.boundary(r['gate_ids'])[0 if side=='prefix' else 1]}
            self.assertEqual(expected,actual)
