"""A published curve requires its observed epochs and selected checkpoint."""

from pathlib import Path
import sys
from unittest import TestCase
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from evaluate_loss import check_trace


class TraceTests(TestCase):
    def test_validation_minimum_binds_to_the_selected_checkpoint(self):
        trace = pd.DataFrame({"epoch": [1, 2, 3], "validation_loss": [3., 1., 2.],
                              "validation_mse": [1., 2., 3.]})
        trained = {"best_epoch": 2, "epochs": 3, "validation_loss": 1., "validation_mse": 2.}
        config = {"max_epochs": 3, "patience": 30}
        self.assertEqual(check_trace(trace, trained, config), 2)
        with self.assertRaisesRegex(ValueError, "selected checkpoint"):
            check_trace(trace, {**trained, "best_epoch": 1}, config)
