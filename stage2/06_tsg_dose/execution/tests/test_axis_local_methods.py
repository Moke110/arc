"""Observe successor dose interfaces on the accepted native artifact."""
import sys
import unittest

import numpy as np
from pathlib import Path

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / 'output/api'))
from tsg_round1.learner import ScientificFailure, serialize_prediction
from tsg_round1.model import configuration, fit_transform


class MethodTests(unittest.TestCase):
    def test_configuration_matches_compact_fixed_protocol(self):
        for location, features, prediction in [('input', 2080, 500), ('latent', 2048, 532)]:
            for vector in ('numeric', 'embedding'):
                with self.subTest(location=location, vector=vector):
                    method = f'{location}-{vector}'
                    transform = fit_transform(method, [0.1, 1., 10.])
                    configuration_ = configuration(method, transform)
                    self.assertEqual(configuration_['dose_vector_dim'], 32)
                    self.assertEqual(configuration_['features_dim'], features)
                    self.assertEqual(configuration_['predictive_head_dim'], prediction)
                    self.assertEqual(configuration_['dose_embedding'],
                                     [1, 32, 32] if vector == 'embedding' else None)
                    self.assertAlmostEqual(transform['mean'], 0.)
                    self.assertAlmostEqual(transform['std'], (2 / 3) ** .5)

    def test_unrepresentable_prediction_is_a_scientific_failure(self):
        for value in (float('nan'), float('inf'), 1e100):
            with self.subTest(value=value), self.assertRaisesRegex(ScientificFailure, 'non-finite'):
                serialize_prediction(np.array([value], dtype=np.float64))
