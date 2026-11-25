"""Verify the formal report's exact native control boundary."""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest import TestCase

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from report_sweep import verify_native_control


class NativeControlTests(TestCase):
    def test_matches_native_columns_and_rejects_a_changed_native_component(self):
        with TemporaryDirectory() as scratch:
            current, previous = (Path(scratch) / name for name in ("fresh", "previous"))
            native_components = (
                "control_mse", "treatment_mse", "delta", "control_kl", "treatment_kl", "prediction_kl"
            )
            fields = ["epoch", "train_loss", "validation_loss", "validation_mse"]
            fields.extend(f"{part}_{component}" for part in ("train", "validation")
                          for component in native_components)
            fresh_fields = [name + "_mse" if name.endswith("_delta") else name for name in fields]
            fresh_fields.extend(("train_delta_l1", "validation_delta_l1"))
            values = list(range(1, len(fields) + 1))
            for run, columns, row in ((current, fresh_fields, values + [70, 90]),
                                     (previous, fields, values)):
                run.mkdir()
                (run / "learning_curve.csv").write_text(
                    ",".join(columns) + "\n" + ",".join(map(str, row)) + "\n"
                )
                (run / "training.json").write_text(json.dumps({"best_epoch": 1, "epochs": 1}))
                for name in ("validation_prediction.npy", "prediction.npy",
                             "predicted_delta.npy", "measured_delta.npy"):
                    (run / name).write_bytes(b"exact fixture array")
            self.assertEqual(verify_native_control(current, previous)["best_epoch"], 1)

            values[4] = 0
            (previous / "learning_curve.csv").write_text(
                ",".join(fields) + "\n" + ",".join(map(str, values)) + "\n"
            )
            with self.assertRaisesRegex(ValueError, "native objective trace"):
                verify_native_control(current, previous)
