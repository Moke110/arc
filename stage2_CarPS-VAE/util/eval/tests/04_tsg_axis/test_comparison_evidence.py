"""Scoring a sealed prediction must use the registered matched plate control."""

from hashlib import sha256
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str((Path(__file__).resolve().parents[4] / "04_tsg_axis/execution/tests" / Path(__file__).name).parents[2] / "output/api"))

from tsg_axis.inputs import AxisInputs
from tsg_axis.comparison_evidence import read_fixed_jaccard


class FixedJaccardTest(unittest.TestCase):
    def test_matched_delta_is_reconstructed_from_fixed_prediction(self):
        control = np.array([[10, 20]], dtype=np.float32)
        delta = np.array([[1, 0]], dtype=np.float32)
        data = AxisInputs(SimpleNamespace(
            manifest_sha256="release",
            manifest={"inputs": {"axis_sha256": "axis"}},
            genes=pd.DataFrame({"gene_index": [0, 1], "symbol": ["A", "B"]}),
            treatments=pd.DataFrame({
                "plate": ["P1"], "drug_parent": ["D1"], "dose_um": [1.0],
                "control_row": [0], "cold_drug_0": ["test"],
            }),
            controls=["P1"],
            read_controls=lambda start, stop: control[start:stop],
            read_treatments=lambda start, stop, layer: delta[start:stop],
        ))
        with tempfile.TemporaryDirectory() as scratch:
            run = Path(scratch)
            data.genes.to_csv(run / "genes.csv", index=False)
            data.condition_table([0]).to_csv(run / "conditions.csv", index=False)
            np.save(run / "prediction.npy", control + delta)
            np.save(run / "measured_delta.npy", delta)
            np.save(run / "predicted_delta.npy", delta)
            config = {
                "training_sha256": "release", "axis_sha256": "axis",
                "partitions": {"test": [0]}, "conditions": {"test": list(data.condition_ids([0]))},
            }

            def manifest():
                return {"fold": 0, "test_conditions": 1, "outputs": {
                    path.name: sha256(path.read_bytes()).hexdigest() for path in run.iterdir()
                }}

            conditions, top200, weighted = read_fixed_jaccard(run, config, manifest(), data)
            self.assertEqual(conditions, data.condition_ids([0]))
            np.testing.assert_array_equal(top200, [1.0])
            np.testing.assert_array_equal(weighted, [1.0])

            np.save(run / "predicted_delta.npy", np.array([[0, 1]], dtype=np.float32))
            with self.assertRaisesRegex(ValueError, "matched-control predicted_delta"):
                read_fixed_jaccard(run, config, manifest(), data)
