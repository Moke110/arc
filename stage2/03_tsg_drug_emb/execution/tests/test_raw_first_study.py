"""Exercise the representation study's raw-first comparison boundary."""
import importlib.util
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
DOSE = STEP.parent / "04_tsg_dose"
sys.path[:0] = [str(STEP / "output/api"), str(DOSE / "output/api")]
from tsg_raw_first import preparation
from tsg_raw_first.inputs import Inputs

spec = importlib.util.spec_from_file_location(
    "paired_fixtures", DOSE / "execution/tests/test_raw_first.py"
)
fixtures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixtures)


class TestRawFirstStudy(unittest.TestCase):
    def test_prepare_score_and_render_pinned_representation(self):
        from stage2_util import metrics, palette, plotting
        from tsg_drug_emb_raw_first import Comparison, evaluate, prepare, render

        with TemporaryDirectory() as temp:
            root = Path(temp)
            representations = fixtures.make_representations(root)
            pair = root / "pair"
            preparation.prepare(pair, view=fixtures.make_views(root / "view"),
                                representations=representations)
            for source in ("tenv", "tahoe_bulk"):
                with self.subTest(source=source):
                    data = Inputs(pair / source)
                    directory, digest = representations[source]
                    receipt = root / f"{source}.json"
                    prepare(data.directory, receipt, representation=directory,
                            representation_sha256=digest)
                    comparison = Comparison(
                        receipt, data.directory, representation=directory, representation_sha256=digest
                    )
                    self.assertEqual(comparison.features("KPGT", [0, 3])[:, 0].tolist(), [-4., -3.])
                    rows = data.rows("cold-drug", 0, "test")
                    predictions = np.array(data.expression[rows])
                    ids = data.obs.index[rows]
                    report = root / f"{source}-report"
                    with patch.object(metrics, "vector_metrics", wraps=metrics.vector_metrics) as called:
                        evaluate(comparison, "KPGT", 0, predictions, ids, report)
                        self.assertEqual(called.call_count, len(rows))
                    summary = pd.read_csv(report / "summary.csv")
                    self.assertEqual(summary.loc[0, "pcc"], 1.)
                    self.assertEqual(summary.loc[0, "parent_equal_pcc"], 1.)
                    self.assertEqual(summary.loc[0, "representation_sha256"], digest)
                    self.assertEqual(summary.loc[0, "source"], source)
                    with self.assertRaisesRegex(ValueError, "prediction identity"):
                        evaluate(comparison, "KPGT", 0, predictions, ids[::-1], root / "bad")
                    csv_only = root / f"{source}.csv"
                    csv_only.write_bytes((report / "summary.csv").read_bytes())
                    png = root / f"{source}.png"
                    with (patch.object(plotting, "read_rows", wraps=plotting.read_rows) as read,
                          patch.object(palette, "color", wraps=palette.color) as color):
                        render(csv_only, png)
                        read.assert_called_once_with(csv_only)
                        self.assertTrue(color.called)
                    self.assertTrue(png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
                    receipt.write_text("{}")
                    with self.assertRaisesRegex(ValueError, "receipt differs"):
                        Comparison(receipt, data.directory, representation=directory,
                                   representation_sha256=digest)
