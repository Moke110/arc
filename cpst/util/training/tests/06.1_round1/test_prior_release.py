from copy import deepcopy
from itertools import product
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst/06_tsg_shape/06.1_round1/execution"))
from evaluate_shape import validate_coordinates, verify_configurations
from prior_shape import load_prior, PRIOR_RECEIPT


class PriorReleaseTest(unittest.TestCase):
    def test_report_rejects_an_unsealed_receipt(self):
        with TemporaryDirectory() as temporary:
            destination = Path(temporary)
            receipt = destination / "received.json"
            altered = json.loads(PRIOR_RECEIPT.read_text())
            altered["fits"] = 30
            receipt.write_text(json.dumps(altered))
            with self.assertRaisesRegex(ValueError, "sealed #347 receipt"):
                load_prior(destination / "matrix", receipt)

    def test_report_requires_each_coordinate_once(self):
        shapes = ("original", "widen_z-110", "widen_z-120", "widen_z-130",
                  "widen_enc-2k", "widen_enc-4k", "widen_enc-6k", "deepen_enc", "deepen_dec", "deepen_both")
        expected = tuple(product(shapes, range(3)))
        validate_coordinates(expected, expected)
        with self.assertRaisesRegex(ValueError, "duplicate shape coordinate"):
            validate_coordinates(expected, (*expected, expected[0]))
        with self.assertRaisesRegex(ValueError, "incomplete shape comparison"):
            validate_coordinates(expected, expected[:-1])

    @unittest.skipUnless("TSG_PRIOR_MATRIX" in os.environ, "requires the materialized #347 matrix")
    def test_report_loads_the_accepted_fit_identities_under_current_names(self):
        paths, fits, provenance = load_prior(Path(os.environ["TSG_PRIOR_MATRIX"]), PRIOR_RECEIPT)
        self.assertEqual(len(paths), 24)
        self.assertEqual(list(fits), list(paths))
        for fold in range(3):
            self.assertEqual(paths["original", fold].parts[-2], "widen_z-100")
            self.assertEqual(paths["widen_enc-6k", fold].parts[-2], "widen_enc")
            self.assertEqual(paths["deepen_both", fold].parts[-2], "deepen-both")
        self.assertEqual(provenance["release_revision"], "9dee5fd90ea6d2bed775655b0205875c8812938b")

    @unittest.skipUnless("TSG_PRIOR_MATRIX" in os.environ, "requires the materialized #347 matrix")
    def test_report_rejects_an_incompatible_optimizer_in_a_mixed_comparison(self):
        matrix = Path(os.environ["TSG_PRIOR_MATRIX"])
        path = matrix / "runs/widen_z-100/fold-0"
        accepted = json.loads((path / "registration.json").read_text())
        candidate = deepcopy(accepted)
        candidate["model"].update(encoder=[2000], learning_rate=0.001)
        candidates = {
            ("original", 0): (path, accepted, None, None, None),
            ("widen_enc-2k", 0): (path, candidate, None, None, None),
        }
        with self.assertRaisesRegex(ValueError, "mismatched protocol or populations"):
            verify_configurations(candidates)
        candidate["model"]["learning_rate"] = 0.0001
        verify_configurations(candidates)

    def test_report_rejects_a_matrix_manifest_that_does_not_match_the_prior_receipt(self):
        with TemporaryDirectory() as temporary:
            matrix = Path(temporary)
            (matrix / "manifest.json").write_text('{"fits": []}')
            with self.assertRaisesRegex(ValueError, "sealed #347 source"):
                load_prior(matrix, PRIOR_RECEIPT)
