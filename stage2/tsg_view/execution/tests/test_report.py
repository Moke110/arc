"""Report-table contract: aggregation identities, provenance and deterministic CSVs.

The fixture writes a complete 45-origin run tree, a frozen-truth tree, diagnostics and a
synthetic accepted ``Mean_cell`` table.  Every value is a literal formula in this file,
so the tests recompute the expected equal-fold mean, fold SD, seed SD and paired interval
without calling the implementation under test.
"""
import hashlib
import io
import json
import math
import shutil
import statistics
import sys
import tempfile
import unittest
from pathlib import Path
from contextlib import redirect_stdout

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import report  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
SOURCE_ORDER = ("tenv", "tahoe_bulk")
VIEW_ORDER = {
    "tenv": ("full", "merge", "solo", "full-pooled-label", "solo-deepest"),
    "tahoe_bulk": ("full", "merge", "solo"),
}
DOSES = {"tenv": (0.2, 1.0, 10.0), "tahoe_bulk": (0.05, 0.5, 5.0)}
SAMPLES = {"tenv": (4, 6, 8), "tahoe_bulk": (5, 7, 9)}
CONDITIONS = {"tenv": (2, 2, 3), "tahoe_bulk": (2, 3, 3)}
VIEW_BASE = {"full": 0.30, "merge": 0.36, "solo": 0.33, "full-pooled-label": 0.40,
             "solo-deepest": 0.38}
METRIC_BASE = {"pcc": 0.20, "spearman": 0.24, "rmse": 0.44, "r2": -0.50}
NOISE_FLOOR = {"tenv": 0.3425617959827319, "tahoe_bulk": 0.07025318359737137}
NOISE_MEAN = {"tenv": 0.318634595773231, "tahoe_bulk": 0.10259702766174636}
CEILING = {"tenv": 0.3802732892389079, "tahoe_bulk": 0.42828477670427056}
CONDITIONS_TOTAL = {"tenv": 1377, "tahoe_bulk": 9741}
PROFILES_TOTAL = {"tenv": 4206, "tahoe_bulk": 22282}


def seeds(view):
    return (42,) if view == "solo-deepest" else (42, 43)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def sample_value(source, view, fold, seed, index, metric):
    """The literal per-profile metric formula the fixture writes."""
    shift = 0.03 if source == "tahoe_bulk" else 0.0
    base = VIEW_BASE[view] + 0.01 * fold + 0.002 * (seed - 42) + shift
    if metric == "pcc":
        return base + 0.005 * index
    if metric == "spearman":
        return 0.9 * base + 0.004 * index
    if metric == "rmse":
        return 0.30 + 0.01 * fold + 0.002 * (seed - 42) + shift + 0.003 * index
    if metric == "r2":
        return METRIC_BASE["r2"] + 0.02 * fold + 0.003 * (seed - 42) + 0.05 * list(
            VIEW_ORDER[source]).index(view) + 0.001 * index
    raise AssertionError(metric)


def condition_value(source, view, fold, seed, metric):
    shift = 0.02 if source == "tahoe_bulk" else 0.0
    base = 0.20 + 0.01 * fold + 0.001 * (seed - 42) + shift
    offsets = {"pcc": 0.0, "spearman": 0.01, "rmse": 0.20, "r2": -0.80}
    return base + offsets[metric] + 0.02 * list(VIEW_ORDER[source]).index(view)


def _summary(samples):
    result = {"n_samples": int(len(samples))}
    for metric in report.METRICS:
        values = samples[metric].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        result[metric] = float(finite.mean()) if len(finite) else None
        result[f"undefined_{metric}_count"] = int(len(values) - len(finite))
    result["constant_truth_count"] = int(samples["constant_truth"].sum())
    result["constant_prediction_count"] = int(samples["constant_prediction"].sum())
    result["nonfinite_count"] = int(samples["nonfinite"].sum())
    return result


def _run(samples, origin, n, *, nan_metrics=()):
    seed = origin["seed"]
    view = origin["view"]
    for index in range(n):
        row = {"profile_id": f"{origin['source']}::{view}::f{origin['fold']}::s{seed}::p{index}",
               "drug_parent": "PARENT", "dose_um": DOSES[origin["source"]][origin["fold"]],
               "cell_line": "CL", "constant_truth": False, "constant_prediction": False,
               "nonfinite": False}
        for metric in report.METRICS:
            value = sample_value(origin["source"], view, origin["fold"], seed, index, metric)
            if metric in nan_metrics and index == 0:
                value = float("nan")
            row[metric] = value
        samples.append(row)
    return pd.DataFrame(samples)


def build_fixture(root, *, run_sources=("tenv", "tahoe_bulk"), nan_origins=()):
    """Write a synthetic run tree, truth tree, diagnostics and Mean_cell provenance."""
    root = Path(root)
    runs = root / "runs"
    truth = root / "truth"
    mean_cell_path = root / "accepted" / "round1_runs.csv"
    receipt_path = root / "accepted" / "round1-baselines-acceptance.json"
    for source in SOURCE_ORDER:
        base = truth / source
        write_json(base / "manifest.json", {"schema": "tsg-view-truth-v1", "source": source,
                                            "scenario": "cold-drug-dose"})
        write_json(base / "diagnostics.json", {
            "schema": "tsg-view-truth-diagnostics-v1", "source": source,
            "scenario": "cold-drug-dose", "label": "heuristic diagnostic",
            "noise_floor": {"median": NOISE_FLOOR[source], "mean": NOISE_MEAN[source],
                            "definition": "fixture"},
            "pcc_ceiling": {"value": CEILING[source], "label": "heuristic diagnostic",
                            "definition": "fixture", "input_conditions": CONDITIONS_TOTAL[source],
                            "input_profiles": PROFILES_TOTAL[source]}})
        for fold in range(3):
            dose = DOSES[source][fold]
            n = SAMPLES[source][fold]
            m = CONDITIONS[source][fold]
            primary = pd.DataFrame({
                "profile_id": [f"{source}::p{fold}-{index}" for index in range(n)],
                "cell_line": "CL", "drug_parent": "PARENT", "dose_um": dose,
                "partition": "test", "fold": fold, "scenario": "cold-drug-dose"})
            primary_path = base / "primary" / f"fold-{fold}" / "population.csv"
            primary_path.parent.mkdir(parents=True, exist_ok=True)
            primary.to_csv(primary_path, index=False, lineterminator="\n")
            counts = [n // m + (1 if index < n % m else 0) for index in range(m)]
            if source == "tahoe_bulk" and fold == 0:
                # like the real Tahoe truth, the replicated-condition subset is smaller
                # than the primary test population
                counts[-1] -= 1
            ids = list(primary["profile_id"])
            members, start = [], 0
            for count in counts:
                members.append(ids[start:start + count])
                start += count
            condition = pd.DataFrame({
                "profile_id": [f"{source}::c{fold}-{index}" for index in range(m)],
                "cell_line": "CL", "drug_parent": "PARENT", "dose_um": dose,
                "members": [";".join(member) for member in members],
                "count": counts, "partition": "test"})
            condition_path = base / "condition_mean" / f"fold-{fold}" / "population.csv"
            condition_path.parent.mkdir(parents=True, exist_ok=True)
            condition.to_csv(condition_path, index=False, lineterminator="\n")
    baseline_rows = []
    for source in SOURCE_ORDER:
        for fold in range(3):
            for scenario in ("cold-drug", "cold-drug-dose"):
                for method in ("Mean_cell", "ECFP4_MLP"):
                    row = {"source": source, "scenario": scenario, "seed": 42,
                           "method": method, "fold": fold, "n_samples": SAMPLES[source][fold]}
                    for metric in report.METRICS:
                        row[metric] = (METRIC_BASE[metric] + 0.01 * fold
                                       + (0.03 if source == "tahoe_bulk" else 0.0)
                                       + (0.5 if method != "Mean_cell" else 0.0))
                    baseline_rows.append(row)
    mean_cell_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(baseline_rows).to_csv(mean_cell_path, index=False, lineterminator="\n")
    write_json(receipt_path, {"schema": "tsg-dose-round1-baselines-report-v1",
                              "status": "passed",
                              "files": {"data/round1_runs.csv": sha256(mean_cell_path)}})
    for origin in report.plan_origins():
        if origin["source"] not in run_sources:
            continue
        key = (origin["source"], origin["view"], origin["fold"], origin["seed"])
        run = runs / origin["source"] / origin["view"] / f"fold-{origin['fold']}" / (
            f"seed-{origin['seed']}")
        run.mkdir(parents=True, exist_ok=True)
        write_json(run / "config.json", {"origin": origin, "mode": "formal",
                                         "truth": {}, "freeze": {}})
        samples = _run([], origin, SAMPLES[origin["source"]][origin["fold"]],
                       nan_metrics=("pcc",) if key in nan_origins else ())
        samples.to_csv(run / "test_per_sample.csv", index=False, lineterminator="\n")
        primary = _summary(samples)
        condition = {"n_samples": CONDITIONS[origin["source"]][origin["fold"]]}
        for metric in report.METRICS:
            condition[metric] = condition_value(origin["source"], origin["view"],
                                               origin["fold"], origin["seed"], metric)
            condition[f"undefined_{metric}_count"] = 0
        metrics = dict(primary)
        metrics.update({"primary": primary, "condition_mean": condition,
                        "n_primary_samples": primary["n_samples"],
                        "n_condition_mean_samples": condition["n_samples"]})
        write_json(run / "test_metrics.json", metrics)
        epochs = [{"epoch": epoch,
                   "train_loss": 1.0 / (epoch + 1) + 0.001 * len(origin["view"]),
                   "validation_mse": 2.0 / (epoch + 2) + 0.002 * origin["fold"]}
                  for epoch in range(1, 5 + origin["fold"] + 1)]
        pd.DataFrame(epochs).to_csv(run / "learning_curve.csv", index=False,
                                    lineterminator="\n")
        write_json(run / "terminal.json", {"status": "complete", "fold": origin["fold"],
                                           "seed": origin["seed"],
                                           "test_samples": primary["n_samples"]})
        files = {name: sha256(run / name) for name in (
            "config.json", "terminal.json", "test_metrics.json", "test_per_sample.csv",
            "learning_curve.csv")}
        write_json(run / "manifest.json", {"schema": "tsg-view-run-v1",
                                           "origin": origin, "outputs": files})
    return {"root": root, "runs": runs, "truth": truth, "mean_cell": mean_cell_path,
            "receipt": receipt_path, "mean_cell_sha": sha256(mean_cell_path)}


def report_args(fixture, out):
    return ["--runs", str(fixture["runs"]), "--truth", str(fixture["truth"]),
            "--report", str(Path(out) / "report"),
            "--figure-data", str(Path(out) / "figure"),
            "--root", str(fixture["root"]),
            "--mean-cell-source", str(fixture["mean_cell"]),
            "--mean-cell-sha256", fixture["mean_cell_sha"],
            "--mean-cell-receipt", str(fixture["receipt"])]


def build(fixture, out):
    result = report.build_report(
        runs_root=fixture["runs"], truth_root=fixture["truth"], repo_root=fixture["root"],
        mean_cell_source=fixture["mean_cell"], mean_cell_sha256=fixture["mean_cell_sha"],
        mean_cell_receipt=fixture["receipt"])
    report.write_report(result, Path(out) / "report", Path(out) / "figure")
    return result


class ReportFixture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temporary = tempfile.TemporaryDirectory(prefix="tsg-view-report-")
        cls.tmp = Path(cls._temporary.name)
        cls.fixture = build_fixture(cls.tmp / "fixture")
        cls.out = cls.tmp / "out"
        cls.result = build(cls.fixture, cls.out)

    @classmethod
    def tearDownClass(cls):
        cls._temporary.cleanup()

    def table(self, name):
        return pd.read_csv(self.out / "report" / f"{name}.csv",
                           float_precision="round_trip")

    def figure(self, name):
        return pd.read_csv(self.out / "figure" / f"{name}.csv",
                           float_precision="round_trip")


class TestAggregationIdentities(ReportFixture):
    def test_primary_rows_are_run_sample_equal_scores(self):
        table = self.table("metrics_primary")
        self.assertEqual(len(table), 45)
        for row in table.to_dict("records"):
            n = SAMPLES[row["source"]][int(row["fold"])]
            self.assertEqual(row["n_samples"], n)
            values = [sample_value(row["source"], row["view"], int(row["fold"]),
                                   int(row["seed"]), index, "pcc") for index in range(n)]
            self.assertAlmostEqual(row["delta_pcc"], float(np.mean(values)), places=12)
            self.assertEqual(row["arm_type"],
                             "solo-deepest sensitivity" if row["view"] == "solo-deepest"
                             else "main")

    def test_equal_fold_mean_uses_unweighted_fold_means(self):
        primary = self.table("metrics_primary")
        folds = self.table("fold_spread")
        primary = primary[primary["truth"] == "primary"]
        folds = folds[folds["truth"] == "primary"]
        self.assertEqual(len(folds), 60)
        for key, block in primary.groupby(["source", "view", "seed"]):
            observed = folds[(folds["source"] == key[0]) & (folds["view"] == key[1])
                             & (folds["seed"] == key[2]) & (folds["metric"] == "delta_pcc")]
            self.assertEqual(len(observed), 1)
            fold_means = [float(block[block["fold"] == fold]["delta_pcc"].iloc[0])
                          for fold in sorted(block["fold"])]
            self.assertAlmostEqual(float(observed["equal_fold_mean"].iloc[0]),
                                   float(np.mean(fold_means)), places=12)
            self.assertAlmostEqual(float(observed["fold_sd"].iloc[0]),
                                   float(statistics.stdev(fold_means)), places=12)
            counted = [int(block[block["fold"] == fold]["n_defined_pcc"].iloc[0])
                       for fold in sorted(block["fold"])]
            pooled = float(np.average(fold_means, weights=counted))
            self.assertAlmostEqual(float(observed["sample_equal_pooled"].iloc[0]), pooled,
                                   places=12)

    def test_sample_equal_pooled_differs_from_equal_fold_headline(self):
        primary = self.table("metrics_primary")
        seeds = self.table("seed_spread")
        seeds = seeds[(seeds["truth"] == "primary") & (seeds["metric"] == "delta_pcc")]
        primary = primary[primary["truth"] == "primary"]
        head = seeds[(seeds["source"] == "tenv") & (seeds["view"] == "merge")]
        block = primary[(primary["source"] == "tenv") & (primary["view"] == "merge")]
        self.assertAlmostEqual(float(head["headline"].iloc[0]),
                               float(block["delta_pcc"].mean()), places=12)
        self.assertNotAlmostEqual(float(head["sample_equal_pooled"].iloc[0]),
                                  float(head["headline"].iloc[0]), places=6)

    def test_seed_spread_reports_two_seeds_separately(self):
        folds = self.table("fold_spread")
        seeds = self.table("seed_spread")
        for view in ("full", "merge", "solo", "full-pooled-label"):
            block = folds[(folds["truth"] == "primary") & (folds["source"] == "tenv")
                          & (folds["view"] == view) & (folds["metric"] == "delta_pcc")]
            per_seed = [float(block[block["seed"] == seed]["equal_fold_mean"].iloc[0])
                        for seed in (42, 43)]
            row = seeds[(seeds["truth"] == "primary") & (seeds["source"] == "tenv")
                        & (seeds["view"] == view) & (seeds["metric"] == "delta_pcc")]
            self.assertEqual(int(row["n_seeds"].iloc[0]), 2)
            self.assertEqual(int(row["n_folds_total"].iloc[0]), 3)
            self.assertAlmostEqual(float(row["headline"].iloc[0]), float(np.mean(per_seed)),
                                   places=12)
            self.assertAlmostEqual(float(row["seed_sd"].iloc[0]),
                                   float(statistics.stdev(per_seed)), places=12)
            self.assertIn("never pooled", row["estimand"].iloc[0])

    def test_solo_deepest_is_tagged_and_single_seed(self):
        seeds = self.table("seed_spread")
        row = seeds[(seeds["source"] == "tenv") & (seeds["view"] == "solo-deepest")
                    & (seeds["truth"] == "primary") & (seeds["metric"] == "delta_pcc")]
        self.assertEqual(int(row["n_seeds"].iloc[0]), 1)
        self.assertTrue(math.isnan(float(row["seed_sd"].iloc[0])))
        self.assertEqual(row["arm_type"].iloc[0], "solo-deepest sensitivity")

    def test_undefined_metric_is_dropped_not_zeroed(self):
        with tempfile.TemporaryDirectory(prefix="tsg-view-report-undefined-") as temporary:
            root = Path(temporary)
            fixture = build_fixture(root / "fixture",
                                    nan_origins={("tenv", "merge", 0, 42)})
            out = root / "out"
            build(fixture, out)
            primary = pd.read_csv(out / "report" / "metrics_primary.csv")
            row = primary[(primary["source"] == "tenv") & (primary["view"] == "merge")
                          & (primary["fold"] == 0) & (primary["seed"] == 42)]
            self.assertEqual(int(row["n_defined_pcc"].iloc[0]),
                             SAMPLES["tenv"][0] - 1)
            values = [sample_value("tenv", "merge", 0, 42, index, "pcc")
                      for index in range(1, SAMPLES["tenv"][0])]
            self.assertAlmostEqual(float(row["delta_pcc"].iloc[0]),
                                   float(np.mean(values)), places=12)
            folds = pd.read_csv(out / "report" / "fold_spread.csv")
            fold_row = folds[(folds["truth"] == "primary") & (folds["source"] == "tenv")
                             & (folds["view"] == "merge") & (folds["seed"] == 42)
                             & (folds["metric"] == "delta_pcc")]
            fold_means = [float(row["delta_pcc"].iloc[0])]
            block = primary[(primary["source"] == "tenv") & (primary["view"] == "merge")
                            & (primary["seed"] == 42) & (primary["fold"] > 0)]
            fold_means.extend(float(value) for value in block["delta_pcc"])
            self.assertAlmostEqual(float(fold_row["equal_fold_mean"].iloc[0]),
                                   float(np.mean(fold_means)), places=12)


class TestTruthVariants(ReportFixture):
    def test_condition_mean_counts_and_caveat(self):
        table = self.table("metrics_condition_mean")
        self.assertEqual(len(table), 45)
        for row in table.to_dict("records"):
            self.assertEqual(int(row["fold"]), int(row["fold"]))
            self.assertEqual(int(row["n_conditions"]),
                             CONDITIONS[row["source"]][int(row["fold"])])
            self.assertEqual(int(row["n_profiles"]),
                             SAMPLES[row["source"]][int(row["fold"])])
            self.assertEqual(int(row["truth_replicated_condition_count"]),
                             CONDITIONS[row["source"]][int(row["fold"])])
            self.assertEqual(row["caveat_flag"], True)
            self.assertIn("not a random subset", row["caveat"])
            self.assertEqual(row["truth"], "condition_mean")

    def test_per_dose_panels_are_keyed_to_the_held_dose(self):
        table = self.table("metrics_per_dose")
        primary = self.table("metrics_primary")
        self.assertEqual(len(table), 45 * 4)
        for source in SOURCE_ORDER:
            block = table[table["source"] == source]
            self.assertEqual(sorted({round(value, 6) for value in block["held_dose_um"]}),
                             sorted(DOSES[source]))
            for fold in range(3):
                fold_block = block[block["fold"] == fold]
                self.assertEqual({round(value, 6) for value in fold_block["held_dose_um"]},
                                 {round(DOSES[source][fold], 6)})
                self.assertEqual({int(value) for value in fold_block["n_samples"]},
                                 {SAMPLES[source][fold]})
        for row in table.to_dict("records"):
            self.assertEqual(row["dose_fold_confounded"], True)
            self.assertIn("confounded by design", row["design_note"])
            match = primary[(primary["source"] == row["source"])
                            & (primary["view"] == row["view"])
                            & (primary["fold"] == row["fold"])
                            & (primary["seed"] == row["seed"])]
            self.assertAlmostEqual(row["value"], float(match[row["metric"]].iloc[0]), places=12)

    def test_diagnostics_beside_scores(self):
        table = self.table("diagnostics_beside_scores")
        seeds = self.table("seed_spread")
        self.assertEqual(len(table), 8)
        for row in table.to_dict("records"):
            self.assertEqual(row["noise_floor"], NOISE_FLOOR[row["source"]])
            self.assertEqual(row["pcc_ceiling"], CEILING[row["source"]])
            self.assertEqual(row["pcc_ceiling_label"], "heuristic diagnostic")
            self.assertIn("population SD", row["noise_floor_estimand"])
            expected = seeds[(seeds["truth"] == "primary")
                             & (seeds["source"] == row["source"])
                             & (seeds["view"] == row["view"])
                             & (seeds["metric"] == "delta_pcc")]
            self.assertAlmostEqual(row["delta_pcc"], float(expected["headline"].iloc[0]),
                                   places=12)
            self.assertAlmostEqual(row["pcc_minus_ceiling"],
                                   row["delta_pcc"] - CEILING[row["source"]], places=12)
            self.assertAlmostEqual(row["rmse_minus_noise_floor"],
                                   row["delta_rmse"] - NOISE_FLOOR[row["source"]], places=12)


class TestPairedAndBaseline(ReportFixture):
    def test_paired_difference_interval(self):
        primary = self.table("metrics_primary")
        table = self.table("paired_vs_full")
        self.assertNotIn("solo-deepest", set(table["view"]))
        self.assertEqual(len(table), (3 + 2) * 4)
        for row in table.to_dict("records"):
            source, view, metric = row["source"], row["view"], row["metric"]
            pairs = []
            for fold in range(3):
                for seed in (42, 43):
                    mine = primary[(primary["source"] == source) & (primary["view"] == view)
                                   & (primary["fold"] == fold) & (primary["seed"] == seed)]
                    full = primary[(primary["source"] == source) & (primary["view"] == "full")
                                   & (primary["fold"] == fold) & (primary["seed"] == seed)]
                    pairs.append(float(mine[metric].iloc[0]) - float(full[metric].iloc[0]))
            matrix = np.asarray(pairs, dtype=float).reshape(2, 3)
            seed_means = matrix.mean(axis=1)
            headline = float(seed_means.mean())
            fold_sd = float(np.mean(matrix.std(axis=1, ddof=1)))
            seed_sd = float(statistics.stdev(seed_means))
            half = 1.96 * math.sqrt((fold_sd / math.sqrt(3)) ** 2
                                    + (seed_sd / math.sqrt(2)) ** 2)
            self.assertEqual(int(row["n_pairs"]), 6)
            self.assertEqual(int(row["missing_pairs"]), 0)
            self.assertAlmostEqual(row["delta_mean"], headline, places=12)
            self.assertAlmostEqual(row["fold_sd"], fold_sd, places=12)
            self.assertAlmostEqual(row["seed_sd"], seed_sd, places=12)
            self.assertAlmostEqual(row["low"], headline - half, places=12)
            self.assertAlmostEqual(row["high"], headline + half, places=12)
            self.assertIn("fold and seed spread are never pooled", row["interval_estimator"])
            self.assertEqual(row["status"], "ok")

    def test_solo_deepest_paired_difference_is_separate(self):
        table = self.table("paired_vs_full_solo_deepest")
        self.assertEqual(set(table["view"]), {"solo-deepest"})
        self.assertEqual(len(table), 4)
        for row in table.to_dict("records"):
            self.assertEqual(int(row["n_pairs"]), 3)
            self.assertEqual(int(row["n_seeds"]), 1)
            self.assertTrue(math.isnan(float(row["seed_sd"])))
            self.assertIn("single-seed sensitivity", row["interval_estimator"])
            self.assertEqual(row["arm_type"], "solo-deepest sensitivity")
            self.assertAlmostEqual(row["low"],
                                   row["delta_mean"] - 1.96 * row["fold_sd"] / math.sqrt(3),
                                   places=12)

    def test_mean_cell_reuses_the_accepted_rows_and_pins_provenance(self):
        table = self.table("baseline_mean_cell")
        self.assertEqual(len(table), 2 * 3 * 4)
        provenance = pd.read_csv(self.fixture["mean_cell"])
        for row in table.to_dict("records"):
            match = provenance[(provenance["source"] == row["source"])
                               & (provenance["scenario"] == "cold-drug-dose")
                               & (provenance["method"] == "Mean_cell")
                               & (provenance["fold"] == row["fold"])]
            self.assertEqual(len(match), 1)
            metric = row["metric"].replace("delta_", "")
            self.assertAlmostEqual(row["value"], float(match[metric].iloc[0]), places=12)
            self.assertEqual(int(row["n_samples"]), SAMPLES[row["source"]][row["fold"]])
            self.assertEqual(row["provenance_sha256"], self.fixture["mean_cell_sha"])
            self.assertEqual(row["acceptance_receipt_sha256"],
                             sha256(self.fixture["receipt"]))
            self.assertEqual(row["truth_rows_match"], True)
            self.assertEqual(row["view"], "all")
            self.assertIn("identical for every view", row["view_scope"])
            self.assertIn("Mean_cell", row["estimand"])

    def test_mean_cell_rejects_a_changed_provenance_file(self):
        with self.assertRaisesRegex(report.ReportError, "provenance digest mismatch"):
            report.build_report(
                runs_root=self.fixture["runs"], truth_root=self.fixture["truth"],
                repo_root=self.fixture["root"], mean_cell_source=self.fixture["mean_cell"],
                mean_cell_sha256="0" * 64, mean_cell_receipt=self.fixture["receipt"])


class TestRunInventory(ReportFixture):
    def test_complete_inventory_and_partial_root(self):
        status = self.table("run_status")
        self.assertEqual(len(status), 45)
        self.assertEqual(set(status["status"]), {"complete"})
        with tempfile.TemporaryDirectory(prefix="tsg-view-report-partial-") as temporary:
            root = Path(temporary)
            fixture = build_fixture(root / "fixture", run_sources=("tenv",))
            out = root / "out"
            build(fixture, out)
            status = pd.read_csv(out / "report" / "run_status.csv")
            self.assertEqual(int((status["status"] == "complete").sum()), 27)
            self.assertEqual(int((status["status"] == "missing").sum()), 18)
            self.assertTrue((status[status["status"] == "missing"]["detail"]
                             == "missing run directory").all())
            primary = pd.read_csv(out / "report" / "metrics_primary.csv")
            self.assertEqual(set(primary["source"]), {"tenv"})
            bars = pd.read_csv(out / "figure" / "metrics_bars.csv")
            self.assertEqual(set(bars["source"]), {"tenv"})
            paired = pd.read_csv(out / "report" / "paired_vs_full.csv")
            tahoe = paired[paired["source"] == "tahoe_bulk"]
            self.assertTrue((tahoe["n_pairs"] == 0).all())
            self.assertTrue((tahoe["status"] == "no-pairs").all())
            self.assertTrue((paired[paired["source"] == "tenv"]["n_pairs"] > 0).all())
            effects = pd.read_csv(out / "figure" / "views_effect.csv")
            self.assertEqual(set(effects["source"]), {"tenv"})

    def test_tampered_run_artifact_fails_closed(self):
        with tempfile.TemporaryDirectory(prefix="tsg-view-report-tamper-") as temporary:
            root = Path(temporary)
            fixture = build_fixture(root / "fixture")
            run = fixture["runs"] / "tenv" / "merge" / "fold-0" / "seed-42"
            samples = pd.read_csv(run / "test_per_sample.csv")
            samples.loc[0, "pcc"] = 0.123456789
            samples.to_csv(run / "test_per_sample.csv", index=False, lineterminator="\n")
            with self.assertRaisesRegex(report.ReportError, "digest mismatch"):
                report.build_report(
                    runs_root=fixture["runs"], truth_root=fixture["truth"],
                    repo_root=fixture["root"], mean_cell_source=fixture["mean_cell"],
                    mean_cell_sha256=fixture["mean_cell_sha"],
                    mean_cell_receipt=fixture["receipt"])

    def test_incomplete_terminal_status_is_reported_not_dropped(self):
        with tempfile.TemporaryDirectory(prefix="tsg-view-report-incomplete-") as temporary:
            root = Path(temporary)
            fixture = build_fixture(root / "fixture")
            run = fixture["runs"] / "tenv" / "solo" / "fold-1" / "seed-43"
            terminal = json.loads((run / "terminal.json").read_text())
            terminal["status"] = "validation-complete"
            (run / "terminal.json").write_text(json.dumps(terminal))
            out = root / "out"
            build(fixture, out)
            status = pd.read_csv(out / "report" / "run_status.csv")
            row = status[(status["source"] == "tenv") & (status["view"] == "solo")
                         & (status["fold"] == 1) & (status["seed"] == 43)]
            self.assertEqual(row["status"].iloc[0], "incomplete")
            self.assertIn("validation-complete", row["detail"].iloc[0])
            primary = pd.read_csv(out / "report" / "metrics_primary.csv")
            self.assertEqual(len(primary), 44)


class TestDeterminism(ReportFixture):
    def test_report_csv_bytes_are_identical_on_rerun(self):
        first, second = self.tmp / "rerun-a", self.tmp / "rerun-b"
        with redirect_stdout(io.StringIO()):
            report.main(report_args(self.fixture, first))
            report.main(report_args(self.fixture, second))
        for name in ("report", "figure"):
            first_files = sorted((first / name).glob("*.csv"))
            second_files = sorted((second / name).glob("*.csv"))
            self.assertEqual([path.name for path in first_files],
                             [path.name for path in second_files])
            for left, right in zip(first_files, second_files):
                self.assertEqual(left.read_bytes(), right.read_bytes(), left.name)
        manifest = json.loads((first / "report" / "report_manifest.json").read_text())
        self.assertEqual(manifest["schema"], "tsg-view-report-v1")
        self.assertEqual(manifest["complete_origins"], 45)
        self.assertEqual(manifest["pending_origins"], 0)
        self.assertIn("metrics_primary.csv", manifest["tables"])
        self.assertIn("figure/metrics_bars.csv", manifest["figure_data"])

    def test_figure_data_is_a_strict_projection(self):
        bars = self.figure("metrics_bars")
        seeds = self.table("seed_spread")
        self.assertEqual(len(bars), 8 * 4)
        for row in bars.to_dict("records"):
            match = seeds[(seeds["truth"] == "primary") & (seeds["source"] == row["source"])
                          & (seeds["view"] == row["view"]) & (seeds["metric"] == row["metric"])]
            self.assertAlmostEqual(row["headline"], float(match["headline"].iloc[0]), places=12)
            if int(row["n_seeds"]) == 1:
                self.assertTrue(math.isnan(float(row["seed_sd"])))
                self.assertTrue(math.isnan(float(match["seed_sd"].iloc[0])))
            else:
                self.assertAlmostEqual(row["seed_sd"], float(match["seed_sd"].iloc[0]),
                                       places=12)
        points = self.figure("metrics_points")
        self.assertEqual(len(points), 45 * 4)
        primary = self.table("metrics_primary")
        for row in points.head(10).to_dict("records"):
            match = primary[(primary["source"] == row["source"]) & (primary["view"] == row["view"])
                            & (primary["fold"] == row["fold"])
                            & (primary["seed"] == row["seed"])]
            self.assertAlmostEqual(row["value"], float(match[row["metric"]].iloc[0]), places=12)


if __name__ == "__main__":
    unittest.main()
