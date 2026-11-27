"""Exercise the current Step 08 CSV-only figure entry point."""

import csv
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chem_figures import replay
from chem_runs import ARMS
from cpst_util.comparison import CURVE_FIELDS, METRIC_FIELDS, METRIC_NAMES, write_paired_metrics


class ComparisonTest(unittest.TestCase):
    def test_canonical_figure_replay_uses_only_plot_csvs(self):
        with TemporaryDirectory() as temporary:
            work = Path(temporary)
            data = work / "data"
            data.mkdir()
            for name, fields, rows in (
                ("learning-curve.csv", CURVE_FIELDS, [(arm, 0, 17, 1, 3, 4) for arm in ARMS]),
                ("metrics.csv", METRIC_FIELDS, [(arm, 0, 17, metric, 0.4)
                                                for arm in ARMS for metric in METRIC_NAMES]),
            ):
                with (data / name).open("w", newline="") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(fields)
                    writer.writerows(rows)
            write_paired_metrics(data, "KPGT")
            output = work / "figures"
            result = replay(data, output, work, comparison=True)
            self.assertEqual(set(result["png_sha256"]), {
                "learning-curve.png", "metrics.png", "metrics-paired.png",
            })
            self.assertEqual(result["independent_csv_only_renders"], 2)
            self.assertEqual(result, replay(data, output, work, verify=True, comparison=True))
            (output / "obsolete.png").touch()
            with self.assertRaisesRegex(ValueError, "published Step 08 figures"):
                replay(data, output, work, verify=True, comparison=True)
