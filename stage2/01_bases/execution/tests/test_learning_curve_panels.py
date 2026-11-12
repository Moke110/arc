"""Learning-curve figure: one shade per (fold, seed) run in both hue families.

ADR 0063: training is the blue family, validation the green family, and one shade
identifies one ``(fold, seed)`` run in both.  The drawing call returns the plotted
values so the PNG can be rebuilt from the recorded CSV alone.
"""
import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import palette
from bases import plots

FIXTURES = Path(__file__).resolve().parent / "fixtures"
FIXTURE_PNG = FIXTURES / "learning_curve_fixture.png"
FIXTURE_CSV = FIXTURES / "learning_curve_fixture.csv"


def canonical(value):
    """Make NaN and tuples comparable across a CSV round-trip."""
    if isinstance(value, dict):
        return {key: canonical(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return tuple(canonical(item) for item in value)
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return value


def synthetic_models():
    """One learnable model with four runs; the last run was never dispatched."""
    return [{
        "name": "PRnet",
        "runs": [
            {"run": "fold_0-seed_17", "epoch": (1.0, 2.0, 3.0),
             "train": (2.0, 1.5, 1.2), "validation": (2.1, 1.6, 1.3)},
            {"run": "fold_0-seed_42", "epoch": (1.0, 2.0, 3.0),
             "train": (2.2, 1.7, 1.4), "validation": (2.3, 1.8, 1.5)},
            {"run": "fold_1-seed_17", "epoch": (1.0, 2.0, 3.0),
             "train": (2.4, 1.9, 1.6), "validation": (2.5, 2.0, 1.7)},
            {"run": "fold_1-seed_42", "epoch": (), "train": (), "validation": (),
             "reason": "not_dispatched"},
        ],
    }]


class LearningCurveValueTests(unittest.TestCase):
    def setUp(self):
        self.values = plots.learning_curve_values(synthetic_models())
        self.panel = self.values["panels"][0]
        self.runs = {run["series"]: run for run in self.panel["runs"]}

    def test_one_panel_per_learnable_model_with_epoch_on_the_horizontal_axis(self):
        self.assertEqual(self.values["schema"], "stage2-tenv-learning-curve-values-v1")
        self.assertEqual(self.values["kind"], "learning_curve")
        self.assertEqual(self.values["layout"], {"rows": 1, "columns": 1, "panels": ("PRnet",)})
        self.assertEqual(self.panel["panel"], "PRnet")
        self.assertEqual(self.panel["x_label"], "Epoch")
        self.assertEqual(self.panel["y_label"], "Loss")
        first = self.runs["fold_0-seed_17"]
        self.assertEqual(first["epoch"], (1.0, 2.0, 3.0))
        self.assertEqual(first["train"], (2.0, 1.5, 1.2))
        self.assertEqual(first["validation"], (2.1, 1.6, 1.3))

    def test_train_blue_and_validation_green_share_one_shade_per_run(self):
        shades = [run["shade"] for run in self.panel["runs"]]
        self.assertEqual(shades, [2, 4, 6, 8])
        for run, expected in zip(self.panel["runs"], shades):
            self.assertEqual(run["train_hue"], "blue")
            self.assertEqual(run["validation_hue"], "green")
            self.assertEqual(run["shade"], expected)
            self.assertEqual(palette.run_shade(self.panel["runs"].index(run)), expected)

    def test_undefined_run_is_recorded_as_undefined_with_a_reason(self):
        missing = self.runs["fold_1-seed_42"]
        self.assertFalse(missing["defined"])
        self.assertEqual(missing["reason"], "not_dispatched")
        self.assertEqual(missing["epoch"], ())
        self.assertTrue(all(run["defined"] for run in self.panel["runs"][:3]))

    def test_duplicate_runs_or_non_monotonic_epochs_are_rejected(self):
        models = synthetic_models()
        models[0]["runs"][1]["run"] = "fold_0-seed_17"
        with self.assertRaisesRegex(ValueError, "unique"):
            plots.learning_curve_values(models)
        models = synthetic_models()
        models[0]["runs"][0]["epoch"] = (1.0, 1.0, 3.0)
        with self.assertRaisesRegex(ValueError, "increase"):
            plots.learning_curve_values(models)
        models = synthetic_models()
        models[0]["runs"][0]["train"] = (2.0, float("nan"), 1.2)
        with self.assertRaisesRegex(ValueError, "finite"):
            plots.learning_curve_values(models)


class LearningCurveDrawingTests(unittest.TestCase):
    def figure(self):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        return plt.subplots(1, 1, figsize=(7.5, 4.8), layout="constrained")

    def test_draw_returns_exactly_the_values_it_drew(self):
        values = plots.learning_curve_values(synthetic_models())
        fig, axes = self.figure()
        try:
            drawn = plots.draw_learning_curve(values, axes)
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)
        self.assertEqual(canonical(drawn), canonical(values))

    def test_each_run_draws_blue_training_and_green_validation_at_its_shade(self):
        values = plots.learning_curve_values(synthetic_models())
        fig, axes = self.figure()
        try:
            plots.draw_learning_curve(values, axes)
            panel = axes
            self.assertEqual(panel.lines[0].get_color(), palette.color("blue", 2))
            self.assertEqual(panel.lines[1].get_color(), palette.color("green", 2))
            self.assertEqual(tuple(panel.lines[0].get_xdata()), (1.0, 2.0, 3.0))
            self.assertEqual(tuple(float(value) for value in panel.lines[0].get_ydata()), (2.0, 1.5, 1.2))
            self.assertEqual(panel.lines[-2].get_color(), palette.color("blue", 6))
            self.assertEqual(panel.lines[-1].get_color(), palette.color("green", 6))
            self.assertTrue(any("undefined" in text.get_text() for text in panel.texts))
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)


class LearningCurveReplayTests(unittest.TestCase):
    def test_recorded_values_round_trip_and_replay_byte_identical(self):
        import matplotlib
        matplotlib.use("Agg")
        values = plots.learning_curve_values(synthetic_models())
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            values_csv = tmp / "values.csv"
            plots.write_plot_values(values, values_csv)
            self.assertEqual(canonical(plots.read_plot_values(values_csv)), canonical(values))
            first, second = tmp / "first.png", tmp / "second.png"
            plots.render_learning_curve(synthetic_models(), first)
            plots.replay_learning_curve(values_csv, second)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_committed_fixture_is_the_synthetic_values_and_replays_byte_identical(self):
        import matplotlib
        matplotlib.use("Agg")
        self.assertTrue(FIXTURE_PNG.exists(), "run the fixture build command to create learning_curve_fixture.png")
        self.assertTrue(FIXTURE_CSV.exists(), "run the fixture build command to create learning_curve_fixture.csv")
        self.assertEqual(canonical(plots.read_plot_values(FIXTURE_CSV)),
                         canonical(plots.learning_curve_values(synthetic_models())))
        with tempfile.TemporaryDirectory() as tmp:
            replay = Path(tmp) / "replay.png"
            plots.replay_learning_curve(FIXTURE_CSV, replay)
            self.assertEqual(replay.read_bytes(), FIXTURE_PNG.read_bytes())


if __name__ == "__main__":
    unittest.main()
