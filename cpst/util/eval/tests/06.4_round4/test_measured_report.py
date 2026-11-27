"""Verify the registered Step 06.4 report with known matched observations."""

from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst/06_tsg_shape/06.4_round4/execution"))

from measured_scores import paired_report
from prior_measured import RECEIPT, ROUND_THREE, verify_release, verify_source_scores


SHAPES = ("XX-XX", "XX-XXΔ", "X-XΔ", "XΔ-XΔ", "XΔ-XΔΔ")
METRICS = ("Δ PCC", "Δ RMSE", "Δ Spearman", "Δ E-distance")


class PairedReportTest(unittest.TestCase):
    def setUp(self):
        self.coordinates = tuple(
            (shape, fold, seed) for shape in SHAPES for fold in (0, 1) for seed in (17, 42))
        values = {
            "XX-XX": (1, 2, 10, 11), "XX-XXΔ": (2, 3, 11, 12), "X-XΔ": (3, 4, 12, 13),
            "XΔ-XΔ": (8, 5, 19, 19), "XΔ-XΔΔ": (4, 6, 15, 18),
        }
        self.scores = pd.DataFrame([
            {"shape": shape, "fold": fold, "seed": seed, "metric": metric,
             "value": values[shape][2 * fold + (seed == 42)],
             "n_conditions": 4, "n_eligible": 4, "n_undefined": 0}
            for shape, fold, seed in reversed(self.coordinates) for metric in METRICS
        ])

    def test_seed_pairs_weight_equally_in_each_fold_and_overall(self):
        _, _, pairs, folds, overall = paired_report(self.scores, self.coordinates, METRICS)
        self.assertEqual((len(pairs), len(folds), len(overall)), (32, 16, 8))
        target = "XΔ-XΔΔ − XX-XXΔ"
        seed_pairs = pairs.loc[pairs.contrast.eq(target) & pairs.metric.eq("Δ Spearman")]
        self.assertEqual(seed_pairs.value.tolist(), [2, 3, 4, 6])
        fold_means = folds.loc[folds.contrast.eq(target) & folds.metric.eq("Δ Spearman")]
        self.assertEqual(fold_means.value.tolist(), [2.5, 5.0])
        summary = overall.loc[overall.contrast.eq(target) & overall.metric.eq("Δ Spearman")]
        self.assertEqual(summary.value.item(), 3.75)
        self.assertAlmostEqual(summary.sample_sd.item(), 1.707825127659933)
        sharing = pairs.loc[pairs.contrast.eq("XΔ-XΔ − XΔ-XΔΔ") & pairs.metric.eq("Δ E-distance")]
        self.assertEqual(sharing.value.tolist(), [4, -1, 4, 1])

    def test_missing_duplicate_and_mismatched_coordinates_fail(self):
        for scores in (
            self.scores.iloc[1:],
            pd.concat([self.scores, self.scores.iloc[[0]]], ignore_index=True),
            self.scores.assign(seed=self.scores.seed.where(self.scores.seed == 17, 43)),
            self.scores.assign(fold=self.scores.fold.where(self.scores.fold == 1, 2)),
        ):
            with self.subTest(rows=len(scores), seeds=set(scores.seed), folds=set(scores.fold)):
                with self.assertRaises(ValueError):
                    paired_report(scores, self.coordinates, METRICS)

    def test_undefined_score_invalidates_only_pairs_using_it(self):
        missing = (self.scores["shape"].eq("XX-XXΔ") & self.scores.seed.eq(17)
                   & self.scores.fold.eq(0) & self.scores.metric.eq("Δ PCC"))
        self.scores.loc[missing, ["value", "n_eligible", "n_undefined"]] = [float("nan"), 0, 4]
        _, _, _, _, summary = paired_report(self.scores, self.coordinates, METRICS)
        pcc = summary.loc[summary.metric.eq("Δ PCC")]
        self.assertEqual(pcc.n_defined.tolist(), [3, 4])



class SourceReportTest(unittest.TestCase):
    def test_source_manifest_must_match_accepted_receipt(self):
        scratch = REPO / "tmp"
        scratch.mkdir(exist_ok=True)
        with TemporaryDirectory(prefix="step064-source-test-", dir=scratch) as directory:
            candidate = Path(directory)
            receipt = candidate / RECEIPT
            receipt.parent.mkdir()
            shutil.copyfile(ROUND_THREE / RECEIPT, receipt)
            altered = candidate / "output/formal/matrix/manifest.json"
            altered.parent.mkdir(parents=True)
            altered.write_text('{"fits": []}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "source artifact differs from the accepted release"):
                verify_release(candidate)

    def test_reused_score_and_eligibility_must_match(self):
        accepted = pd.DataFrame([{"shape": "XX-XX", "fold": 0, "seed": 17, "metric": "Δ PCC",
                                 "run_sha256": "registered", "value": .2, "n_conditions": 4,
                                 "n_eligible": 4, "n_undefined": 0}])
        bandwidth = pd.DataFrame([{"fold": 0, "bandwidth": .3}])
        verify_source_scores(accepted.copy(), bandwidth.copy(), accepted, bandwidth)
        with self.assertRaisesRegex(ValueError, "rescored source predictions or eligibility differ"):
            verify_source_scores(accepted.assign(n_eligible=3, n_undefined=1), bandwidth, accepted, bandwidth)
        with self.assertRaisesRegex(ValueError, "rescored source predictions or eligibility differ"):
            verify_source_scores(accepted.assign(run_sha256="other"), bandwidth, accepted, bandwidth)


if __name__ == "__main__":
    unittest.main()
