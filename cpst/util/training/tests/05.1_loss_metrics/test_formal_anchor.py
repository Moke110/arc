"""The native control must replay the accepted Step 04 prediction arrays."""

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest import TestCase

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst/05_tsg_loss/05.1_loss_metrics/execution"))

from report_formal import verify_native_control


class NativeControlTest(TestCase):
    def test_a_one_byte_prediction_change_fails_native_control_acceptance(self):
        with TemporaryDirectory() as scratch:
            root = Path(scratch)
            current, previous = root / "candidate", root / "step04"
            for run in (current, previous):
                run.mkdir()
                (run / "learning_curve.csv").write_text(
                    "epoch,train_loss,validation_loss,validation_mse\n1,1.0,1.0,0.1\n")
                (run / "training.json").write_text('{"best_epoch": 1, "epochs": 1}')
                for name in ("validation_prediction.npy", "prediction.npy",
                             "predicted_delta.npy", "measured_delta.npy"):
                    (run / name).write_bytes(b"exact")
            self.assertIn("prediction.npy", verify_native_control(current, previous)["arrays_sha256"])
            (current / "prediction.npy").write_bytes(b"exacx")
            with self.assertRaisesRegex(ValueError, "prediction.npy"):
                verify_native_control(current, previous)
