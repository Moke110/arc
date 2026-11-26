"""Verify a completed matrix can be reopened without changing its audit."""

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matrix_composition as matrix
from composition_roster import coordinates, fit_coordinates


class MatrixRerunTest(unittest.TestCase):
    def test_completed_matrix_is_read_only(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            output, audit = root / "matrix", root / "audit"
            audit.mkdir()
            samples = audit / "resource-samples.jsonl"
            samples.write_text("previous acceptance measurement\n")
            expected, fits = {}, []
            for key, path in matrix.run_paths(output, True).items():
                shape, fold, seed = key
                config = expected[key] = {"shape": shape, "fold": fold, "seed": seed}
                arguments = {**config, "axis": matrix.AXIS, "device": "cuda", "smoke": True}
                registration = matrix.axis_run.write_json(path / "registration.json",
                                                          {**config, "arguments": arguments, "nonce": repr(key)})
                run = matrix.axis_run.write_json(path / "manifest.json", {"status": "PASS"})
                fits.append({**config, "path": path.relative_to(output).as_posix(),
                             "registration_sha256": matrix.axis_run.digest(registration),
                             "run_manifest_sha256": matrix.axis_run.digest(run)})
            evidence = {"coordinates": [list(key) for key in coordinates(True)],
                        "new_fits": [list(key) for key in fit_coordinates(True)]}
            preflight = matrix.axis_run.write_json(output / "preflight.json", evidence)
            receipt = matrix.axis_run.write_json(output / "manifest.json", {"schema": matrix.SCHEMA, "fits": fits,
                                                  "preflight_sha256": matrix.axis_run.digest(preflight)})
            receipt_sha = matrix.axis_run.digest(receipt)
            with patch.object(matrix, "preflight", return_value=(expected, evidence)):
                result = matrix.execute(output, audit, root, root, root, True)
            self.assertEqual(result, {"result": "PASS", "fits": 24, "manifest_sha256": receipt_sha})
            self.assertEqual(list(audit.iterdir()), [samples])
            self.assertEqual(samples.read_text(), "previous acceptance measurement\n")
            self.assertEqual(matrix.axis_run.digest(receipt), receipt_sha)


if __name__ == "__main__":
    unittest.main()
