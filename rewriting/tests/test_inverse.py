import random
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

from prefix_classifier.model import Gate
from prefix_classifier.inverse import angle_key,dagger_gate,inverse_relation,cancel_inverses,schemas,replay_inverse_certificate
from prefix_classifier.qcanopy import sequence_to_dag
from prefix_classifier import load_library
from prefix_classifier.pair import prove_pair

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'rules/verified'))
import exact


def dag(n,gs,phase=0):
    return sequence_to_dag(n,[dict(name=g.op,qubits=list(g.qubits),params=list(g.params)) for g in gs],phase)


def symbolic_exp(angle):
    key=angle_key(angle)
    if key is None:raise ValueError(angle)
    exps=[0,0,0,0]
    for name,value in key:
        if name=='':
            if value:raise ValueError('non-cyclotomic numeric constant')
        elif name=='pi':
            power=4*value
            if power.denominator!=1:raise ValueError(angle)
            exps[0]+=int(power)
        else:
            if value.denominator!=1:raise ValueError(angle)
            exps[1+('theta','phi','lam').index(name)]+=int(value)
    return exact.polynomial([(tuple(exps),1)])


def matrix(g):
    if g.op=='csxdg':
        # Definition of controlled SX†: identity on control 0, SX† on control 1.
        b=exact.gate_matrix('sxdg')
        return tuple(tuple((exact.ONE if r==c else exact.ZERO) if not c&1 else
            (b[r>>1][c>>1] if r&1 else exact.ZERO) for c in range(4)) for r in range(4))
    return exact.gate_matrix(g.op,tuple(g.params))


def product(a,b):
    return tuple(tuple(sum_poly(exact.mul(a[r][k],b[k][c]) for k in range(len(a)))
        for c in range(len(a))) for r in range(len(a)))


def sum_poly(values):
    out=exact.ZERO
    for x in values:out=exact.add(out,x)
    return out


class InverseSchemaTests(unittest.TestCase):
    def test_every_schema_has_exact_symbolic_matrix_proof(self):
        # Full formal matrices, not sampled angles or floating tolerances.
        rule=schemas();gs=[]
        for op,(_,n) in rule['fixed'].items():gs.append(Gate('a',op,tuple(range(n)),()))
        for op,n in rule['negate_angle'].items():gs.append(Gate('a',op,tuple(range(n)),('theta',)))
        gs.append(Gate('a','u',(0,),('theta','phi','lam')))
        exact.gate_matrix.cache_clear()
        try:
            with patch.object(exact,'exp_half',side_effect=symbolic_exp):
                for g in gs:
                    inv=dagger_gate(g)
                    self.assertIsNotNone(inverse_relation(g,inv))
                    a,b=matrix(g),matrix(inv)
                    self.assertEqual(product(b,a),exact.identity(len(a)),g.op)
                    self.assertEqual(product(a,b),exact.identity(len(a)),g.op)
        finally:exact.gate_matrix.cache_clear()

    def test_exact_numeric_and_affine_parameters_no_snapping(self):
        a=Gate('a','rx',(0,),(.123456789,))
        self.assertIsNotNone(inverse_relation(a,Gate('b','rx',(0,),(-.123456789,))))
        self.assertIsNone(inverse_relation(a,Gate('b','rx',(0,),(-.123456789+2**-50,))))
        self.assertEqual(angle_key('-(theta+2*phi-pi/2)'),angle_key('-theta-2*phi+pi/2'))
        self.assertIsNotNone(inverse_relation(Gate('a','rz',(0,),('pi/2',)),Gate('b','rz',(0,),('-0.5*pi',))))
        self.assertIsNone(inverse_relation(Gate('a','rx',(0,),('pi',)),Gate('b','rx',(0,),(-3.141592653589793,))))
        for x in ('theta*phi','sin(theta)','__import__("os")','nan'):
            self.assertIsNone(angle_key(x))

    def test_wrong_ports_u_parameter_order_and_unknown_gates(self):
        self.assertIsNone(inverse_relation(Gate('a','cx',(0,1)),Gate('b','cx',(1,0))))
        g=Gate('a','u',(0,),(.1,.2,.3))
        self.assertIsNotNone(inverse_relation(g,Gate('b','u',(0,),(-.1,-.3,-.2))))
        self.assertIsNone(inverse_relation(g,Gate('b','u',(0,),(-.1,-.2,-.3))))
        self.assertIsNone(dagger_gate(Gate('a','unknown',(0,))))
        self.assertIsNone(dagger_gate(Gate('a','unitary_1q',(0,),(1,0,0,0,0,0,1,0))))

    def test_groups_cancel_inside_out_with_disjoint_interleaving(self):
        gs=[Gate('a','h',(0,)),Gate('b','cx',(0,1)),Gate('c','ry',(1,),('pi/4',))]
        noise=Gate('z','z',(2,))
        d=dag(3,gs+[noise]+[dagger_gate(g) for g in reversed(gs)],phase=.125)
        before=d.to_dict();r=cancel_inverses(d)
        self.assertEqual(len(r.witnesses),3)
        self.assertEqual([g.op for g in r.core_dag.nodes],['z'])
        self.assertEqual(r.core_dag.global_phase,.125)
        self.assertEqual(before,d.to_dict())
        self.assertTrue(replay_inverse_certificate(d,r.to_dict()))

    def test_interactions_and_wrong_group_order_block_cancellation(self):
        gs=[Gate('a','h',(0,)),Gate('b','x',(0,))]
        wrong=dag(1,gs+[dagger_gate(g) for g in gs])
        self.assertEqual(len(cancel_inverses(wrong).core_dag.nodes),4)
        d=dag(2,[Gate('a','x',(0,)),Gate('b','cx',(0,1)),Gate('c','x',(0,))])
        self.assertEqual(len(cancel_inverses(d).core_dag.nodes),3)

    def test_random_gate_groups_preserve_exact_full_unitary(self):
        rng=random.Random(759);exact.gate_matrix.cache_clear();exact.unitary.cache_clear()
        try:
            with patch.object(exact,'exp_half',side_effect=symbolic_exp):
                for _ in range(16):
                    gs=[]
                    for i in range(8):
                        op=rng.choice(['h','x','s','t','sx','cx','cz','rx','rz','u'])
                        qs=(0,1) if op in ('cx','cz') else (rng.randrange(2),)
                        ps=('pi/4','-pi/2','pi') if op=='u' else (('pi/4',) if op in ('rx','rz') else ())
                        gs.append(Gate(str(i),op,qs,ps))
                    d=dag(2,gs+[dagger_gate(g) for g in reversed(gs)])
                    r=cancel_inverses(d)
                    self.assertFalse(r.core_dag.nodes)
                    self.assertTrue(replay_inverse_certificate(d,r.to_dict()))
                    steps=tuple((d.by_id[g].op,d.by_id[g].qubits,d.by_id[g].params) for g in d.order)
                    self.assertEqual(exact.unitary(2,steps),exact.unitary(2,()))
        finally:exact.gate_matrix.cache_clear();exact.unitary.cache_clear()

    def test_long_group_before_pair_search_cap_and_explicit_opt_out(self):
        lib=load_library(ROOT/'rules/verified/classes.phase.json')
        gs=[Gate(str(i),'rx',(0,),(i+.125,)) for i in range(100)]
        d=dag(1,gs+[dagger_gate(g) for g in reversed(gs)])
        r=prove_pair(dag(1,[]),d,lib,max_residual_gates=16)
        self.assertEqual(r['verdict'],'equivalent_up_to_global_phase')
        self.assertTrue(any(x['kind']=='inverse_identity' for x in r['proof']))
        self.assertEqual(prove_pair(dag(1,[]),d,lib,max_residual_gates=16,enable_inverse=False)['verdict'],'inconclusive')

    def test_near_inverse_must_not_prove_identity(self):
        lib=load_library(ROOT/'rules/verified/classes.phase.json')
        d=dag(1,[Gate('a','rx',(0,),(.125,)),Gate('b','rx',(0,),(-.125+2**-40,))])
        self.assertFalse(cancel_inverses(d).witnesses)
        self.assertEqual(prove_pair(dag(1,[]),d,lib)['verdict'],'inconclusive')

    def test_certificate_replay_rejects_nonadjacency_and_duplicate_consumption(self):
        import copy
        good=dag(1,[Gate('a','x',(0,)),Gate('b','x',(0,))])
        cert=cancel_inverses(good).to_dict(include_core=False)
        self.assertTrue(replay_inverse_certificate(good,cert))
        bad=copy.deepcopy(cert);bad['witnesses']+=copy.deepcopy(bad['witnesses'])
        self.assertFalse(replay_inverse_certificate(good,bad))
        separated=dag(1,[Gate('a','x',(0,)),Gate('b','z',(0,)),Gate('c','x',(0,))])
        from prefix_classifier.model import digest
        bad=copy.deepcopy(cert);bad['input_fingerprint']=digest(separated.to_dict())
        bad['witnesses'][0]['lhs_gate_ids']=['g0','g2']
        self.assertFalse(replay_inverse_certificate(separated,bad))
