"""Round 1 score semantics and the validation-to-test gate."""

import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd

from test_gene_axes import build
from test_landmark import make_fixture, run

SCRIPT = Path(__file__).resolve().parents[1] / "script"
sys.path.insert(0, str(SCRIPT))
from round1_fit import run_queue
from round_contract import role_result, sha256, verify_committed_promotion
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
        rows = [{"model": "XPert_base", "auprc_full": .4}, {"model": "XPert_G2000", "auprc_full": .4},
                {"model": "XPert_G5000", "auprc_full": .6}, {"model": "XPert_MLP512", "auprc_full": .7}]
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


def commit(repo, *files):
    subprocess.run(["git", "-C", str(repo), "add", *files], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=fixture@example.invalid",
                    "-c", "user.name=Fixture", "commit", "-qm", "seal"], check=True)


def write_labels(round_root, corpus, role, promotion=None):
    labels = round_root / "labels" / role
    labels.mkdir(parents=True)
    full = corpus / MODEL_VIEWS["XPert_G18164"]
    landmark = corpus / MODEL_VIEWS["XPert_base"]
    genes = pd.read_csv(full / "genes.csv")
    supplement = pd.read_csv(landmark / "genes.csv").query("symbol == 'HDGFL3'")
    pd.concat([genes, supplement])[["hgnc_id", "symbol"]].to_csv(labels / "genes.csv", index=False)
    groups = pd.read_csv(full / f"tenv/fold_0/{role}.csv").group_id.drop_duplicates().tolist()
    groups.reverse()
    pd.DataFrame({"group_id": groups, "status": "eligible", "reason": ""}).to_csv(
        labels / "groups.csv", index=False
    )
    (labels / "contrast_samples.csv").write_text("group_id\n")
    padj = np.ones((len(groups), len(genes) + 1))
    padj[:, -1] = .001  # HDGFL3 alone is significant.
    for name, data in (("padj", padj), ("pvalue", padj), ("log2_fold_change", np.full_like(padj, -1.))):
        np.save(labels / f"{name}.npy", data)
    manifest = {"schema": "step09-significance-v1", "role": role, "fold": 0, "fdr": .05,
                "gene_axes_manifest_sha256": sha256(corpus / "gene_axes/manifest.json"),
                "predictor_manifest_sha256": sha256(full / "manifest.json"),
                "evaluator_manifest_sha256": sha256(full.parent / "evaluator/manifest.json"),
                "landmark_predictor_manifest_sha256": sha256(landmark / "manifest.json"),
                "native_sha256": json.loads((corpus / "gene_axes/manifest.json").read_text())["source_hashes"]["tenv.h5ad"],
                "files": {f.name: sha256(f) for f in labels.iterdir()}}
    if promotion:
        manifest["promotion_sha256"] = sha256(promotion)
    (labels / "manifest.json").write_text(json.dumps(manifest))


class RoundWorkflowTest(unittest.TestCase):
    def test_cli_seals_before_test_and_replays_both_figures(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            native, tenv = make_fixture(root)
            axes, landmark = build(root, native, tenv)
            corpus = root / "corpus"
            corpus.mkdir()
            axes.rename(corpus / "gene_axes")
            landmark.rename(corpus / "tahoe_tenv_bulk_l1000")
            subprocess.run(["git", "init", "-q", "-b", "dev", str(root)], check=True)
            policy = root / "policy.json"
            policy.write_text(json.dumps({
                "kind": "qualification", "device": "cpu", "kernel": "math", "dtype": "float32",
                "microbatch": 1, "effective_batch": 1,
                "phases": {"tahoe": {"updates": 1, "validation_interval": 1, "learning_rate": .0001},
                           "tenv": {"updates": 1, "validation_interval": 1, "learning_rate": .0001}},
            }))
            plan = root / "plan.json"
            plan.write_text(json.dumps({
                "schema": "step09-round1-plan-v1", "fold": 0, "seed": 17,
                "policy_sha256": sha256(policy), "references": ["no_change", "delta_mean"],
                "models": [{"name": model, "view": view,
                            "predictor_sha256": sha256(corpus / view / "manifest.json")}
                           for model, view in MODEL_VIEWS.items()],
                "gene_axes_manifest_sha256": sha256(corpus / "gene_axes/manifest.json"),
                "shared_gene_count": 1, "continuous_metrics_provisional": True,
            }))
            commit(root, "plan.json", "policy.json")
            round_root = root / "round1"
            round_root.mkdir()
            work = root / "work"
            model = "XPert_base"
            predictor = corpus / MODEL_VIEWS[model]
            fit = work / "fits" / model
            run("script/xpert.py", "fit", "--predictor", predictor, "--policy", policy, "--fold", 0,
                "--seed", 17, "--output", fit)
            first = work / "resources" / model / "fit-1"
            second = work / "resources" / model / "fit-2"
            first.mkdir(parents=True)
            second.mkdir()
            (first / "time.log").write_text("Elapsed (wall clock) time (h:mm:ss or m:ss): 0:12.00\n")
            (second / "time.log").write_text("Elapsed (wall clock) time (h:mm:ss or m:ss): 0:30.00\n")
            (second / "wall-clock.tsv").write_text("100\t132\n")
            prediction = round_root / "predictions/validation" / model
            run("script/xpert.py", "predict", "--predictor", predictor, "--checkpoint", fit / "tenv.pt",
                "--role", "validation", "--output", prediction)
            write_labels(round_root, corpus, "validation")
            promotion = root / "promotion.json"
            arguments = ("--round-root", round_root, "--work-root", work, "--corpus", corpus, "--plan", plan,
                         "--promotion", promotion)
            run("script/evaluate.py", "seal", *arguments)
            valid = pd.read_csv(round_root / "validation_scores/metrics.csv").set_index("model")
            self.assertTrue(np.isnan(valid.at[model, "auprc_full"]))
            self.assertGreaterEqual(valid.at[model, "auprc_g"], 0.)
            self.assertAlmostEqual(valid.at[model, "fit_attempt_wall_minutes"], 44 / 60)
            self.assertEqual(valid.at[model, "timed_fit_attempts"], 2)
            rejected = run("script/evaluate.py", "test", *arguments, check=False)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("committed dev bytes", rejected.stderr)
            commit(root, "promotion.json")
            write_labels(round_root, corpus, "test", promotion)
            prediction = round_root / "predictions/test" / model
            run("script/xpert.py", "predict", "--predictor", predictor, "--checkpoint", fit / "tenv.pt",
                "--role", "test", "--output", prediction)
            run("script/evaluate.py", "test", *arguments)
            run("script/evaluate.py", "replay", *arguments)
            original_seal = sha256(promotion)
            original_validation = sha256(round_root / "validation_scores/manifest.json")
            run("script/evaluate.py", "revise", *arguments)
            run("script/evaluate.py", "replay-revision", *arguments)
            self.assertEqual(sha256(promotion), original_seal)
            self.assertEqual(sha256(round_root / "legacy_unsigned/validation_scores/manifest.json"),
                             original_validation)
            results = pd.read_csv(round_root / "test_scores/metrics.csv").set_index("model")
            self.assertEqual(results.index.tolist(), ["no_change", "delta_mean", *MODEL_VIEWS])
            self.assertEqual(results.at["XPert_G18164", "status"], "incomplete")
            self.assertEqual(results.at["delta_mean", "training_minutes"], 0.)
            for name in ("xpert_opti_r1_metrics.png", "xpert_opti_r1_PR.png"):
                with (round_root / name).open("rb") as handle:
                    self.assertEqual(handle.read(8), b"\x89PNG\r\n\x1a\n")
