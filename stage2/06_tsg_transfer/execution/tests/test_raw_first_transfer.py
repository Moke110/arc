"""Exercise the active transfer preparation, evaluation and CSV-only renderer."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
DOSE = STEP.parent / "04_tsg_dose"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RawFirstTransferTest(unittest.TestCase):
    def test_transfer_uses_source_only_tahoe_folds_and_frozen_tenv_predictions(self):
        cli = load_module("transfer_raw_first", STEP / "execution/raw_first.py")
        fixtures = load_module("dose_raw_first_fixture", DOSE / "execution/tests/test_raw_first.py")
        from stage2_util import metrics, palette, plotting
        from tsg_raw_first import preparation
        from tsg_raw_first.inputs import Inputs

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pair = root / "pair"
            preparation.prepare(pair, view=fixtures.make_views(root / "view"),
                                representations=fixtures.make_representations(root))
            prepared = root / "prepared"
            cli.main(["prepare", "--pair", str(pair), "--candidate", str(prepared)])
            tahoe = pd.read_csv(prepared / "tahoe/splits.csv")
            self.assertEqual(set(tahoe.fold), {0, 1, 2})
            self.assertTrue(tahoe.profile_id.str.startswith("tahoe_bulk-").all())
            self.assertFalse(tahoe.loc[tahoe.partition == "validation"].empty)

            data = Inputs(pair / "tenv")
            rows = data.rows("cold-drug-dose", 0, "test")
            predictions = root / "predictions.npy"
            np.save(predictions, data.expression[rows])
            identities = root / "identities.csv"
            pd.DataFrame({"profile_id": data.obs.index[rows]}).to_csv(identities, index=False)
            output = root / "report"
            with patch.object(metrics, "vector_metrics", wraps=metrics.vector_metrics) as metric:
                cli.main(["evaluate", "--pair", str(pair), "--prepared", str(prepared), "--fold", "0",
                          "--path", "transfer", "--prediction", str(predictions),
                          "--prediction-ids", str(identities), "--output", str(output)])
                self.assertEqual(metric.call_count, len(rows))
                np.testing.assert_array_equal(metric.call_args.args[0][:3], [1., 2., 5.])
            scores = pd.read_csv(output / "scores.csv")
            self.assertEqual(scores.profile_id.tolist(), data.obs.index[rows].tolist())
            np.testing.assert_allclose(scores.pcc, 1.)

            comparison = root / "comparison"
            paths = []
            for fold in range(3):
                indices = data.rows("cold-drug-dose", fold, "test")
                np.save(predictions, data.expression[indices])
                pd.DataFrame({"profile_id": data.obs.index[indices]}).to_csv(identities, index=False)
                for arm in ("tenv-only", "tahoe-only", "transfer"):
                    assessment = root / f"{arm}-{fold}"
                    cli.main(["evaluate", "--pair", str(pair), "--prepared", str(prepared),
                              "--fold", str(fold), "--path", arm, "--prediction", str(predictions),
                              "--prediction-ids", str(identities), "--output", str(assessment)])
                    paths.append(str(assessment))
            cli.main(["compare", "--output", str(comparison), *paths])
            differences = pd.read_csv(comparison / "transfer_differences.csv")
            self.assertEqual(differences.fold.tolist(), ["0", "1", "2", "mean"])
            np.testing.assert_allclose(differences.transfer_minus_tenv_only, 0.)

            png = root / "out/figure.png"
            with (patch.object(plotting, "read_rows", wraps=plotting.read_rows) as read,
                  patch.object(plotting, "save_figure", wraps=plotting.save_figure) as save,
                  patch.object(palette, "color", wraps=palette.color) as color):
                cli.main(["render", "--data", str(comparison), "--output", str(png)])
                self.assertEqual(read.call_count, 2)
                save.assert_called_once()
                self.assertGreater(color.call_count, 1)
            self.assertEqual(png.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
