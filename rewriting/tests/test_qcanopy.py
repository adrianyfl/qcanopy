import copy
import math
import unittest

from prefix_classifier.model import InputError
from prefix_classifier.qcanopy import parse_qcanopy_ir
from prefix_classifier.swaps import remove_swaps


class QCanopyImportTests(unittest.TestCase):
    def data(self):
        return dict(schema_version=1, n=3, global_phase=.25, gates=[
            dict(name='h', qubits=[0], params=[]),
            dict(name='swap', qubits=[0, 2], params=[]),
            dict(name='cx', qubits=[2, 1], params=[]),
            dict(name='rx', qubits=[1], params=[math.pi / 2])])

    def test_order_phase_parameters_and_mapping(self):
        d = parse_qcanopy_ir(self.data())
        self.assertEqual(d.prev['g2', 2], 'g1')
        self.assertEqual(d.by_id['g2'].qubits, (2, 1))
        self.assertEqual(d.by_id['g3'].params, (math.pi / 2,))
        self.assertEqual(d.global_phase, .25)
        s = remove_swaps(d)
        self.assertEqual(s.output_wire_to_core, (2, 1, 0))
        self.assertEqual(s.core_dag.by_id['g2'].qubits, (0, 1))

    def test_rejects_malformed_and_unsupported(self):
        cases = []
        for key, value in [('n', 33), ('n', True), ('schema_version', True),
                           ('global_phase', float('nan'))]:
            d = self.data(); d[key] = value; cases.append(d)
        for gate in [dict(name='cx', qubits=[0, 0], params=[]),
                     dict(name='h', qubits=[3], params=[]),
                     dict(name='rx', qubits=[0], params=['pi/2']),
                     dict(name='unitary_1q', qubits=[0], params=[]),
                     dict(name='h', qubits=[0], params=[1])]:
            d = self.data(); d['gates'] = [gate]; cases.append(d)
        for d in cases:
            with self.assertRaises(InputError):
                parse_qcanopy_ir(copy.deepcopy(d))
