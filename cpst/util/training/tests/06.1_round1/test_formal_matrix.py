import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst/06_tsg_shape/06.1_round1/execution"))
from matrix_shape import check_registrations, run_paths


class FormalMatrixTest(unittest.TestCase):
    def test_training_requires_every_distinct_preregistered_coordinate(self):
        with TemporaryDirectory() as temporary:
            paths = run_paths(Path(temporary))
            self.assertEqual(list(paths), [(shape, fold) for shape in
                                          ("widen_enc-2k", "widen_enc-4k") for fold in range(3)])
            expected = {}
            for index, ((shape, fold), path) in enumerate(paths.items()):
                identity = {"shape": shape, "fold_identity": fold}
                args = {"shape": shape, "fold": fold, "axis": "denoise-5k",
                        "device": "cuda", "smoke": False}
                record = {**identity, "arguments": args, "nonce": str(index)}
                expected[shape, fold] = identity
                path.mkdir(parents=True)
                (path / "registration.json").write_text(json.dumps(record))
            self.assertEqual(list(check_registrations(paths, expected)), list(paths))
            record["nonce"] = "0"
            (path / "registration.json").write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, "registration differs"):
                check_registrations(paths, expected)
