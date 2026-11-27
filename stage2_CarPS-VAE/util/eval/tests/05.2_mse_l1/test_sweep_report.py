"""Verify the formal report's exact native control boundary."""

from contextlib import ExitStack, contextmanager
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

sys.path.insert(0, str((Path(__file__).resolve().parents[4] / "05_tsg_loss/05.2_mse_l1/execution/tests" / Path(__file__).name).parents[1]))

import report_sweep
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

class PublicationTests(TestCase):
    @contextmanager
    def report_fixture(self, verify=lambda *args: None):
        with TemporaryDirectory() as scratch, ExitStack() as patches:
            root = Path(scratch)
            args = SimpleNamespace(
                **{name: root / name for name in (
                    "matrix", "matrix_audit", "report", "activity", "publication", "scratch",
                    "training", "bandwidth_training", "representations", "source", "previous_step"
                )},
                command="build", code_revision="test revision",
            )
            for directory in (args.matrix, args.matrix_audit, args.report, args.scratch):
                directory.mkdir()
            for directory, names in (
                (args.matrix, ("manifest.json", "preflight.json")),
                (args.matrix_audit, ("resource-samples.jsonl",)),
                (args.report, ("manifest.json", *report_sweep.PROJECTIONS)),
            ):
                for name in names:
                    (directory / name).touch()

            def fixture_render(data, output, activity, attempt):
                output.mkdir(exist_ok=True)
                for name in report_sweep.FIGURES:
                    (output / name).write_text(name)

            for mock in (
                patch.object(report_sweep, "matrix_inputs", return_value=({}, {})),
                patch.object(report_sweep, "command", side_effect=lambda *unused: verify(args)),
                patch.object(report_sweep, "verify_evaluation", return_value={"outputs": {}}),
                patch.object(report_sweep, "verify_predecessor", return_value={}),
                patch.object(report_sweep, "audit_resources", return_value={}),
                patch.object(report_sweep, "render", side_effect=fixture_render),
            ):
                patches.enter_context(mock)
            yield args

    def test_preserves_destination_created_during_verification(self):
        with self.report_fixture(lambda args: args.publication.mkdir()) as args:
            with self.assertRaisesRegex(ValueError, "unused output directory"):
                report_sweep.publish(args)
            self.assertEqual(list(args.publication.iterdir()), [])

    def test_binds_optional_coordinator_samples_using_monotonic_gaps(self):
        with self.report_fixture() as args:
            samples = args.matrix_audit / "coordinator-resource-samples.jsonl"
            samples.write_text(
                '{"seconds_since_start": 0, "unix_seconds": 120, "wsl_available_kib": 500, "gpu_used_mib": 4}\n'
                '{"seconds_since_start": 10.125, "unix_seconds": 100, "wsl_available_kib": 450, "gpu_used_mib": 5}\n'
            )
            report_sweep.publish(args)
            receipt = json.loads((args.publication / "log/333-formal.json").read_text())
            published_samples = "log/333-coordinator-resources.jsonl"
            self.assertEqual((args.publication / published_samples).read_bytes(), samples.read_bytes())
            self.assertIn(published_samples, receipt["published_files_sha256"])
            coordinator = receipt["resources"]["coordinator"]
            self.assertEqual(coordinator["samples"], 2)
            self.assertEqual(coordinator["maximum_observed_sample_gap_seconds"], 10.125)
            self.assertEqual(coordinator["minimum_wsl_available_kib"], 450)
            self.assertEqual(coordinator["maximum_gpu_used_mib"], 5)
