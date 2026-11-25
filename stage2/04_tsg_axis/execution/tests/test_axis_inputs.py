"""The plate-matched release model and fixed-prediction boundaries."""
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / "output/api"))

from tsg_axis.inputs import AxisInputs, score_fixed_prediction


class TrainingReleaseFixture:
    def __init__(self):
        self.genes = pd.DataFrame({"gene_index": [7, 12], "symbol": ["A", "B"]})
        self.controls = pd.DataFrame({"plate": ["p1", "p2"]})
        self.treatments = pd.DataFrame({
            "plate": ["p1", "p2", "p1", "p2"], "drug_parent": ["a", "b", "c", "d"],
            "dose_um": [1., 2., 3., 4.], "control_row": [0, 1, 0, 1],
            "cold_drug_0": ["train", "valid", "test", "test"],
        })
        self.plate = np.asarray([[2, 4], [4, 8]], dtype=np.float32)
        self.expression = np.asarray([[4, 8], [6, 7], [8, 11], [10, 7]], dtype=np.float32)

    def control(self):
        return np.asarray([3, 6], dtype=np.float32)

    def read_controls(self, start, stop):
        return self.plate[start:stop]

    def read_treatments(self, start, stop, layer="treatment"):
        expression = self.expression[start:stop]
        return expression if layer == "treatment" else (
            expression - self.plate[self.treatments.control_row.iloc[start:stop].to_numpy()]
        )


class AxisInputsTest(unittest.TestCase):
    def test_model_partitions_and_delta_arrays_share_matched_plate_controls(self):
        data = AxisInputs(TrainingReleaseFixture())
        train, valid, test = (data.rows(0, phase) for phase in ("train", "valid", "test"))
        np.testing.assert_array_equal(data.model_input("train").control[train], [[2, 4]])
        np.testing.assert_array_equal(data.model_input("valid").control[valid], [[4, 8]])
        np.testing.assert_array_equal(data.model_input("test").control[test], [[2, 4], [4, 8]])
        ids = data.condition_ids(test)
        prediction = np.asarray([[7, 11], [9, 7]], dtype=np.float32)
        result = score_fixed_prediction(data, 0, prediction, ids)
        np.testing.assert_array_equal(result.predicted_delta, [[5, 7], [5, -1]])
        np.testing.assert_array_equal(result.measured_delta, [[6, 7], [6, -1]])
        with self.assertRaisesRegex(ValueError, "condition identities"):
            score_fixed_prediction(data, 0, prediction, ids[::-1])
        with self.assertRaisesRegex(ValueError, "ordered members"):
            score_fixed_prediction(data, 0, prediction, ids, rows=test[::-1])
        with self.assertRaisesRegex(ValueError, "ordered members"):
            score_fixed_prediction(data, 0, prediction[:1], data.condition_ids(train), rows=train)

    def test_prediction_rows_require_exact_integer_positions(self):
        release = TrainingReleaseFixture()
        release.treatments.loc[1, "cold_drug_0"] = "test"
        data = AxisInputs(release)
        ids = data.condition_ids([2, 3])
        prediction = np.ones((2, 2), dtype=np.float32)
        for rows in ([2, 3], [2, np.uint64(3)], np.asarray([2, 3], dtype=np.uint64)):
            with self.subTest(rows=rows):
                self.assertEqual(score_fixed_prediction(data, 0, prediction, ids, rows=rows).condition_ids, ids)
        for rows in ([2.9, 3.9], ["2", "3"], [np.timedelta64(2, "ns"), 3],
                     np.asarray([3, 2], dtype=np.uint64), np.asarray([2, np.iinfo(np.uint64).max], dtype=np.uint64)):
            with self.subTest(rows=rows), self.assertRaisesRegex(ValueError, "ordered members"):
                score_fixed_prediction(data, 0, prediction, ids, rows=rows)
        for boolean in (True, np.bool_(True)):
            with self.subTest(boolean=boolean), self.assertRaisesRegex(ValueError, "ordered members"):
                score_fixed_prediction(data, 0, prediction, data.condition_ids([1, 2]), rows=[boolean, 2])

    def test_selected_gene_width_keeps_native_latent_and_kpgt_dimensions(self):
        import os
        import torch

        dose_api = STEP.parent / "07_tsg_dose/output/api"
        sys.path.insert(0, str(dose_api))
        from tsg_kpgt.model import adapter, fit_transform

        source = Path(os.environ.get("TRANSIGEN_TEST_SOURCE", STEP.parents[1] / "pert_model/TranSiGen/source"))
        if not (source / "src/model.py").is_file():
            self.skipTest("TranSiGen native source is not materialized")
        transform = fit_transform("latent-numeric", [0.2, 1., 10.])
        native = adapter(source, "latent-numeric", transform)
        model = native.build(42, "cpu", n_genes=3000).eval()
        with torch.no_grad():
            outputs = model(torch.zeros((2, 3000)), torch.zeros((2, 2305)))
        self.assertEqual(tuple(outputs[3].shape), (2, 3000))
        self.assertEqual(tuple(outputs[4].shape), (2, 100))
