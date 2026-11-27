"""Exercise the formal sweep preregistration gate on a bounded fixture."""

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest import TestCase

sys.path.insert(0, str((Path(__file__).resolve().parents[4] / "05_tsg_loss/05.2_mse_l1/execution/tests" / Path(__file__).name).parents[1]))

from matrix_sweep import check_registrations, run_paths
from run_sweep import axis_run


class SweepPreregistrationTests(TestCase):
    def test_requires_the_full_ordered_cuda_roster(self):
        with TemporaryDirectory() as scratch:
            paths = run_paths(Path(scratch))
            self.assertEqual(len(paths), 27)
            self.assertEqual(set(paths), {(str(alpha), fold)
                                         for alpha in (0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
                                         for fold in range(3)})
            expected = {key: {"conditions": {"test": ["first", "second"]}, "seed": 42}
                        for key in paths}
            for index, ((alpha, fold), path) in enumerate(paths.items()):
                path.mkdir(parents=True)
                args = dict(alpha=alpha, fold=fold, axis="denoise-5k", device="cuda", smoke=False)
                axis_run.write_json(path / "registration.json",
                                    {**expected[alpha, fold], "arguments": args, "nonce": str(index)})
            self.assertEqual(len(check_registrations(paths, expected)), 27)

            record = paths["2", 2] / "registration.json"
            data = axis_run.read_json(record)
            data["conditions"]["test"].reverse()
            axis_run.write_json(record, data)
            with self.assertRaisesRegex(ValueError, "registration"):
                check_registrations(paths, expected)
