"""Check preregistration completeness and reuse at the formal matrix boundary."""

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest import TestCase

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from matrix import check_registrations, run_paths
from run_loss import axis_run


class FormalRegistrationTests(TestCase):
    def test_rejects_changed_condition_order_before_fitting(self):
        with TemporaryDirectory() as scratch:
            paths = run_paths(Path(scratch))
            self.assertEqual(len(paths), 36)
            self.assertEqual({key[2] for key in paths}, {0, 1, 2})
            self.assertEqual({key[1] for key in paths}, {"1/4", "1/3", "1/2"})
            self.assertEqual({key[0] for key in paths}, {"mse", "cos", "pcc", "l1"})

            expected = {key: {"conditions": {"test": ["first", "second"]}, "seed": 42}
                        for key in paths}
            for index, (key, path) in enumerate(paths.items()):
                path.mkdir(parents=True)
                args = dict(zip(("loss", "alpha", "fold"), key, strict=True), axis="denoise-5k",
                            device="cuda", smoke=False)
                axis_run.write_json(path / "registration.json",
                                    {**expected[key], "arguments": args, "nonce": str(index)})
            self.assertEqual(len(check_registrations(paths, expected)), 36)

            key = ("l1", "1/2", 2)
            record = paths[key] / "registration.json"
            registration = axis_run.read_json(record)
            registration["conditions"]["test"].reverse()
            axis_run.write_json(record, registration)
            with self.assertRaisesRegex(ValueError, "registration"):
                check_registrations(paths, expected)
