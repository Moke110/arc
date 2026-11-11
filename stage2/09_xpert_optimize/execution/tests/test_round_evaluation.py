"""Round 1 score semantics and the validation-to-test gate."""

import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "script"
sys.path.insert(0, str(SCRIPT))
from round1_fit import run_queue
from round_contract import promotion_digest, role_result, sha256, verify_committed_promotion
from round_inputs import MODEL_VIEWS, RoundInputs
from round_metrics import average_precision, continuous, pooled_pr, paired_promotions, project_delta


class RoundEvaluationTest(unittest.TestCase):
    def test_nonfixture_gate_requires_both_roles_terminal(self):
        with TemporaryDirectory() as directory:
            work = Path(directory)
            statuses = work / "status"
            statuses.mkdir()
            inputs = RoundInputs.__new__(RoundInputs)
            inputs.work, inputs.fixture = work, False
            for model in MODEL_VIEWS:
                (statuses / f"{model}.json").write_text(json.dumps({
                    "model": model, "status": "complete", "validation_replay_sha256": "a" * 64,
                }))
            pending = statuses / "XPert_base.json"
            pending.write_text(json.dumps({"model": "XPert_base", "status": "active", "stage": "fit"}))
            with self.assertRaisesRegex(ValueError, "validation requires terminal model statuses: XPert_base"):
                inputs.statuses("validation")
            pending.write_text(json.dumps({"model": "XPert_base", "status": "failed", "stage": "fit", "exit_code": 137}))
            self.assertEqual(len(inputs.statuses("validation")), 9)
            with self.assertRaisesRegex(ValueError, "test requires terminal model statuses: XPert_G2000"):
                inputs.statuses("test")

    def test_default_test_queue_skips_validation_failure_and_records_later_proof_failure(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            work = root / "work"
            statuses = work / "status"
            statuses.mkdir(parents=True)
            failed = {"model": "XPert_base", "status": "failed", "stage": "fit", "exit_code": 137}
            (statuses / "XPert_base.json").write_text(json.dumps(failed))
            later = statuses / "XPert_G2000.json"
            later.write_text(json.dumps({"model": "XPert_G2000", "status": "complete", "attempts": {},
                                         "validation_replay_sha256": "a" * 64}))
            plan = {"models": [{"name": name, "view": name, "predictor_sha256": "b" * 64}
                               for name in ("XPert_base", "XPert_G2000")]}
            self.assertEqual(run_queue(root / "output", work, root / "corpus", plan,
                                       root / "policy.json", "test"), ["XPert_G2000"])
            self.assertEqual(json.loads((statuses / "XPert_base.json").read_text()), failed)
            result = json.loads(later.read_text())
            self.assertEqual(role_result(result, "validation")[0], "complete")
            self.assertEqual(role_result(result, "test")[0], "failed")
            self.assertEqual(result["test_failure"]["stage"], "verify-input-test")

    def test_ties_and_uncovered_positives_keep_fixed_pooled_population(self):
        covered_labels = np.array([0, -1, 0])
        tied = np.array([1., -1., 0.])
        full = average_precision(tied, covered_labels, positives=2)
        model = average_precision(tied, covered_labels, positives=1)
        self.assertEqual((full.ap, model.ap), (.25, .5))
        self.assertEqual((full.recall[-1], model.recall[-1]), (.5, 1.))
        curve = pooled_pr([full, model])
        self.assertTrue(curve.contributing_groups.eq(2).all())
        self.assertAlmostEqual(curve.recall.iloc[-1], 2 / 3)
        self.assertEqual(curve.precision.iloc[-1], .5)
        self.assertEqual(pooled_pr([full]).recall.iloc[-1], .5)

    def test_hdgfl3_and_zero_positives(self):
        self.assertIsNone(average_precision(np.array([0., 0.]), np.zeros(2, dtype=bool), 0))
        self.assertIsNone(average_precision(np.array([0.]), np.array([False]), 0))
        self.assertEqual(average_precision(np.array([0., 2.]), np.array([False, True]), 1).ap, 1.)
        no_hits = average_precision(np.array([2.]), np.array([False]), 1)
        self.assertEqual((no_hits.ap, no_hits.recall[-1]), (0., 0.))
        self.assertEqual(pooled_pr([no_hits]).recall.tolist(), [0.])

    def test_wrong_direction_and_zero_prediction_cannot_hit(self):
        curve = average_precision(np.array([4., -3., 0., 1.]),
                                  np.array([-1, -1, 1, 0]), 3)
        self.assertAlmostEqual(curve.ap, 1 / 6)
        self.assertAlmostEqual(curve.recall[-1], 1 / 3)
        self.assertEqual(curve.precision[-1], 1 / 3)
        silent = average_precision(np.zeros(3), np.ones(3), 3)
        self.assertEqual(silent.ap, 0.)
        self.assertEqual(pooled_pr([silent]).recall.iloc[-1], 0.)
        pooled = pooled_pr([curve, silent])
        self.assertAlmostEqual(pooled.recall.iloc[-1], 1 / 6)
        self.assertTrue(pooled.contributing_groups.eq(2).all())

    def test_shared_projection_uses_log_scale_and_independent_control(self):
        native_control = np.log1p(np.array([[10., 90., 30.], [0., 0., 0.]]))
        prediction = np.array([[0., 0., -100.], [-1., -1., -1.]])
        shared_control = np.log1p(np.array([[8000., 2000.], [0., 0.]]))
        observed = project_delta(prediction, native_control, [0, 2], shared_control)
        np.testing.assert_allclose(observed, [[.22311855, -7.60140233], [0, 0]], rtol=1e-6)

    def test_undefined_correlation_r_squared_and_missing_control(self):
        stats = continuous(np.array([2., 2.]), np.array([2., 2.]))
        self.assertEqual(stats["rmse"], 0.)
        self.assertTrue(all(stats[name] is None for name in ("pcc", "spearman", "r2")))
        rows = [{"model": "XPert_base", "AUPRC": .4}, {"model": "XPert_G2000", "AUPRC": .4},
                {"model": "XPert_G5000", "AUPRC": .6}, {"model": "XPert_MLP512", "AUPRC": .7}]
        decisions = paired_promotions(rows)
        self.assertEqual([item["decision"] for item in decisions], ["not_promoted", "promoted", "unavailable"])

    def test_promotion_gate_checks_committed_dev_bytes(self):
        with TemporaryDirectory() as directory:
            repo = Path(directory)
            subprocess.run(["git", "init", "-q", "-b", "dev", str(repo)], check=True)
            file = repo / "promotion.json"
            record = {"schema": "step09-promotion-v1", "fold": 0, "seed": 17,
                      "plan_sha256": "a" * 64, "validation_manifest_sha256": "b" * 64,
                      "validation_labels_sha256": "c" * 64, "validation_metrics_sha256": "d" * 64,
                      "comparisons": [{}] * 8}
            file.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, "committed dev bytes"):
                verify_committed_promotion(file)
            commit(repo, "promotion.json")
            self.assertEqual(verify_committed_promotion(file), record)
            file.write_text(json.dumps({**record, "seed": 42}))
            with self.assertRaisesRegex(ValueError, "committed dev bytes"):
                verify_committed_promotion(file)

    def test_historical_seal_remains_in_git_without_restoring_score_values(self):
        import hashlib
        with TemporaryDirectory() as directory:
            repo = Path(directory)
            subprocess.run(['git', 'init', '-q', '-b', 'dev', str(repo)], check=True)
            path = repo / 'promotion.json'
            record = {'schema': 'step09-promotion-v1', 'fold': 0, 'seed': 17,
                      **{key: 'a' * 64 for key in ('plan_sha256', 'validation_manifest_sha256',
                                                   'validation_labels_sha256', 'validation_metrics_sha256')},
                      'comparisons': [{'model': f'm{i}', 'control': 'base', 'decision': 'promoted',
                                       'model_auprc_full': .2} for i in range(8)]}
            blob = json.dumps(record).encode()
            path.write_bytes(blob)
            commit(repo, 'promotion.json')
            original = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()
            subprocess.run(['git', '-C', str(repo), 'rm', '-q', 'promotion.json'], check=True)
            pointer = repo / 'promotion-provenance.json'
            metadata = {'schema': 'step09-historical-promotion-v1', 'path': 'promotion.json',
                        'git_commit': original, 'sha256': hashlib.sha256(blob).hexdigest()}
            pointer.write_text(json.dumps(metadata))
            commit(repo, 'promotion-provenance.json')
            self.assertEqual(promotion_digest(path), metadata['sha256'])
            restored = verify_committed_promotion(path)
            self.assertEqual(len(restored['comparisons']), 8)
            self.assertFalse(any('model_auprc_full' in item for item in restored['comparisons']))
            self.assertFalse(path.exists())
            pointer.write_text(json.dumps({**metadata, 'sha256': '0' * 64}))
            commit(repo, 'promotion-provenance.json')
            with self.assertRaisesRegex(ValueError, 'digest'):
                verify_committed_promotion(path)


def commit(repo, *files):
    subprocess.run(["git", "-C", str(repo), "add", *files], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=fixture@example.invalid",
                    "-c", "user.name=Fixture", "commit", "-qm", "seal"], check=True)
