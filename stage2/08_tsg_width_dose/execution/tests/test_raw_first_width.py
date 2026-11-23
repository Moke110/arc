from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


STEP = Path(__file__).resolve().parents[2]
DOSE = STEP.parent / "07_tsg_dose"
CLI = load_module("width_raw_first", STEP / "execution/raw_first.py")
FIXTURE = load_module("dose_raw_first_fixture", DOSE / "execution/tests/test_raw_first.py")


class WidthDoseRawFirstTest(unittest.TestCase):
    def test_prepare_evaluate_and_render_pinned_width_coordinate(self):
        from tsg_raw_first.inputs import KPGTInputs
        from stage2_util import metrics, palette, plotting

        with TemporaryDirectory() as temp:
            root = Path(temp)
            views = FIXTURE.make_views(root / "view")
            representations, digest = FIXTURE.make_representations(root)["tenv"]
            inputs = root / "width-inputs"
            CLI.main(["prepare", "--view", str(views), "--representations", str(representations),
                      "--representation-sha256", digest, "--output", str(inputs)])
            data = KPGTInputs(inputs)
            self.assertEqual(data.config["view_generation"]["view_sha256"],
                             FIXTURE.sha256(views / "manifest.json"))
            self.assertEqual(data.config["training_seed"], 42)
            self.assertEqual(data.fingerprints([0, 3])[:, 0].tolist(), [-4., -3.])

            rows = data.rows("cold-drug-dose", 0, "test")
            prediction = root / "prediction.npy"
            ids = root / "ids.csv"
            report = root / "report"
            np.save(prediction, data.expression[rows])
            pd.DataFrame({"profile_id": data.obs.index[rows]}).to_csv(ids, index=False)
            args = ["evaluate", "--inputs", str(inputs), "--width", "1.5x", "--method", "FiLM",
                    "--fold", "0", "--prediction", str(prediction), "--prediction-ids", str(ids),
                    "--output", str(report)]
            with patch.object(metrics, "vector_metrics", wraps=metrics.vector_metrics) as score:
                CLI.main(args)
                self.assertEqual(score.call_count, len(rows))
            summary = pd.read_csv(report / "summary.csv").iloc[0]
            self.assertEqual((summary.width, summary.method, summary.seed), ("1.5x", "FiLM", 42))
            self.assertEqual(summary.rmse, 0.)
            self.assertEqual(summary.r2, 1.)
            self.assertEqual(summary.pair_sha256, data.pair_sha256)
            self.assertEqual(summary.axis_sha256, data.config["view_generation"]["axis_sha256"])

            summary_path = root / "plot-only.csv"
            summary_path.write_bytes((report / "summary.csv").read_bytes())
            png = root / "metrics.png"
            with (patch.object(plotting, "read_rows", wraps=plotting.read_rows) as read,
                  patch.object(plotting, "save_figure", wraps=plotting.save_figure) as save,
                  patch.object(palette, "color", wraps=palette.color) as color):
                CLI.main(["render", "--summary", str(summary_path), "--output", str(png)])
                read.assert_called_once_with(summary_path)
                save.assert_called_once()
                self.assertGreater(color.call_count, 0)
            self.assertTrue(png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))

            pd.DataFrame({"profile_id": ["other-profile"] * len(rows)}).to_csv(ids, index=False)
            with self.assertRaisesRegex(ValueError, "prediction identity"):
                CLI.main([*args[:-1], str(root / "invalid-report")])
            self.assertFalse((root / "invalid-report").exists())
