"""Public formal matrix contracts around coordinate selection and release reads."""

from pathlib import Path
import json
import shutil
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axis_formal import fresh_coordinates


class FormalRoster(unittest.TestCase):
    def test_unavailable_step08_ledgers_refit_all_nine_references(self):
        from axis_reference import reference_plan

        with TemporaryDirectory() as folder:
            missing = Path(folder) / "missing-step08"
            acceptance = missing / "log/364-formal.json"
            acceptance.parent.mkdir(parents=True)
            for contents in (None, "invalid json", "[]"):
                with self.subTest(contents=contents):
                    if contents is not None:
                        acceptance.write_text(contents)
                    plan = reference_plan(missing, missing, missing, missing)
                    self.assertEqual(len(fresh_coordinates(plan)), 72)
                    self.assertEqual(len(plan["decisions"]), 9)
                    self.assertTrue(all(row["decision"] == "refit" and row["reason"] for row in plan["decisions"]))

    def test_failed_reference_is_fitted_while_successful_replays_are_reused(self):
        decisions = [
            {"coordinate": ["denoise-5k", fold, seed], "decision": "reuse",
             "validation_replay": "exact", "test_replay": "exact"}
            for fold in range(3) for seed in (17, 42, 2026)
        ]
        decisions[3] = {"coordinate": ["denoise-5k", 1, 17], "decision": "refit",
                        "reason": "checkpoint replay differs"}
        remaining = fresh_coordinates({"decisions": decisions})
        self.assertEqual(len(remaining), 64)
        self.assertIn(("denoise-5k", 1, 17), remaining)
        self.assertNotIn(("denoise-5k", 1, 42), remaining)
        with self.assertRaises(ValueError):
            fresh_coordinates({"decisions": decisions[:-1] + decisions[:1]})


class FormalGate(unittest.TestCase):
    def test_gate_rejects_failed_status_and_changed_preflight(self):
        from axis_matrix import verify_bounded_gate
        from axis_protocol import STEP

        with TemporaryDirectory() as folder:
            log = Path(folder)
            for path in (STEP / "log").glob("367-*"):
                shutil.copyfile(path, log / path.name)
            gate = log / "367-bounded.json"
            baseline = json.loads(gate.read_text())
            observed, _ = verify_bounded_gate(log)
            self.assertEqual(len(observed["formal_roster"]), 72)
            gate.write_text(json.dumps({**baseline, "status": "FAIL"}))
            with self.assertRaisesRegex(ValueError, "bounded gate"):
                verify_bounded_gate(log)
            gate.write_text(json.dumps(baseline))
            (log / "367-preflight.json").write_text("{}")
            with self.assertRaisesRegex(ValueError, "bounded gate"):
                verify_bounded_gate(log)


class ScientificAudit(unittest.TestCase):
    def test_plot_audit_rejects_divergent_coverage_or_checkpoint_trace(self):
        import pandas as pd
        from audit_axis_report import audit_plotted_csvs
        from axis_report import figure_data

        with TemporaryDirectory() as folder:
            release = Path(folder)
            run_path = release / "fit"
            run_path.mkdir()
            key = ("hvg-2k", 0, 17)
            identity = dict(zip(("axis", "fold", "seed"), key, strict=True))
            score = {**identity, "scope": "disease", "metric": "direction-correct coverage",
                     "overlap": 1258, "value": 25., "n_eligible": 1, "n_undefined": 0,
                     "paired_difference": -2.}
            summary = {**{name: score[name] for name in ("axis", "scope", "metric", "overlap",
                                                        "n_eligible", "n_undefined")},
                       "mean": 25., "sample_sd": 0.25, "paired_mean": -2., "paired_sd": 0.5}
            objectives = {"epoch": 1, "train_loss": 2.5, "validation_loss": 3.5}
            pd.DataFrame([objectives]).to_csv(run_path / "learning_curve.csv", index=False)
            tables = {"fold-scores.csv": pd.DataFrame([score]), "summary.csv": pd.DataFrame([summary]),
                      "epoch-objectives.csv": pd.DataFrame([{**identity, **objectives}])}
            for name, table in tables.items():
                table.to_csv(release / name, index=False)
            data = release / "plots"
            figure_data(tables, data)
            runs = {key: (run_path,)}
            audit_plotted_csvs(data, release, runs)

            coverage = data / "disease-correct-coverage.csv"
            original = coverage.read_text()
            coverage.write_text(original.replace("25.0", "26.0"))
            with self.assertRaisesRegex(ValueError, "disease-correct-coverage"):
                audit_plotted_csvs(data, release, runs)
            coverage.write_text(original)

            trace = run_path / "learning_curve.csv"
            trace.write_text(trace.read_text().replace("3.5", "3.75"))
            with self.assertRaisesRegex(ValueError, "fit trace"):
                audit_plotted_csvs(data, release, runs)

    def test_audit_rejects_weighted_coverage_and_interchanged_pairs(self):
        from audit_axis_report import audit_aggregates

        scores = {
            ("hvg-2k", 0, 17, "disease", "direction-correct coverage"):
                dict(value="25", reference_value="17", paired_difference="8", overlap="1258",
                     n_conditions="1", n_eligible="1", n_undefined="0"),
            ("hvg-2k", 1, 17, "disease", "direction-correct coverage"):
                dict(value="16", reference_value="18", paired_difference="-2", overlap="1258",
                     n_conditions="5", n_eligible="5", n_undefined="0"),
            ("denoise-5k", 0, 17, "disease", "direction-correct coverage"):
                dict(value="17", reference_value="17", paired_difference="0", overlap="1513",
                     n_conditions="1", n_eligible="1", n_undefined="0"),
            ("denoise-5k", 1, 17, "disease", "direction-correct coverage"):
                dict(value="18", reference_value="18", paired_difference="0", overlap="1513",
                     n_conditions="5", n_eligible="5", n_undefined="0"),
        }
        summary = {
            ("hvg-2k", "disease", "direction-correct coverage"):
                dict(mean="20.5", sample_sd="6.363961030678928", n_defined="2", paired_mean="3",
                     paired_sd="7.0710678118654755", n_paired="2", n_coordinates="2",
                     n_eligible="6", n_undefined="0", overlap="1258"),
            ("denoise-5k", "disease", "direction-correct coverage"):
                dict(mean="17.5", sample_sd="0.7071067811865476", n_defined="2", paired_mean="0",
                     paired_sd="0", n_paired="2", n_coordinates="2",
                     n_eligible="6", n_undefined="0", overlap="1513"),
        }
        audit_aggregates(scores, summary)
        group = summary["hvg-2k", "disease", "direction-correct coverage"]
        group["paired_mean"] = "-0.3333333333333333"
        with self.assertRaises(ValueError):
            audit_aggregates(scores, summary)
        group["paired_mean"] = "3"
        scores["hvg-2k", 0, 17, "disease", "direction-correct coverage"]["reference_value"] = "18"
        with self.assertRaises(ValueError):
            audit_aggregates(scores, summary)
