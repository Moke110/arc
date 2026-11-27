"""Exercise the standalone bounded-similarity renderer with only plot CSVs."""

import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from unittest import TestCase

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
LOSS = REPO / "cpst/05_tsg_loss"
PCC_FIGURE = LOSS / "05.4_mse-l1_pcc/figure"


class SimilarityFigureTests(TestCase):
    def test_renders_the_pcc_study_without_importing_the_training_runtime(self):
        with TemporaryDirectory() as scratch:
            root = Path(scratch)
            (root / "torch.py").write_text("raise RuntimeError('training runtime is unavailable')\n")
            data = root / "data"
            data.mkdir()
            for name in ("delta-metrics.csv", "epoch-objectives.csv", "fold-scores.csv"):
                (data / name).write_bytes((PCC_FIGURE / "data" / name).read_bytes())

            output = root / "figures"
            result = subprocess.run(
                [sys.executable, str(PCC_FIGURE.parent / "execution/render.py"), "--data", str(data),
                 "--output", str(output)],
                env={**os.environ, "PYTHONPATH": str(root)},
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual({path.name for path in output.iterdir()}, {
                "delta-metrics.png", "delta-metrics-improvement.png", "learning-curves.png",
            })
