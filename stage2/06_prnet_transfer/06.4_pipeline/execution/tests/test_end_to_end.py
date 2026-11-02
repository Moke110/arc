"""Fast CPU end-to-end smoke on the synthetic fixture (the ticket's CPU acceptance).

Exercises the real pipeline surface — the 18,164-gene axis, both stages, both
initialisations, strict checkpoint reload, target-free prediction, the frozen
direction-auprc@200 evaluator, model cards — plus the two integrity properties
that are easy to fake: poisoning test labels must not move a training run, and a
perfect prediction must score 1.0 while the trained 2-epoch model scores near
chance on this synthetic data.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch

import _support
from ptlib import checkpoint as checkpoint_mod
from ptlib import data as data_mod
from ptlib import evaluate as evaluate_mod
from ptlib import h5io
from ptlib import modelcard, predict as predict_mod, train as train_mod
from ptlib import model as model_mod

FIXTURE = None
PRETRAIN = None
RANDOM = None
WARM = None
TMP = None


def _pretrain_options(out_dir, seed=5):
    # Hermetic: never resolve the developer's derived cache under output/cache/
    # (the real 4.87 GB memmap lives there and is bound to the frozen source hash).
    cache = Path(out_dir).parent / "cache" / "tahoe_log1p_cp10k_v1.npy"
    return train_mod.TrainOptions(
        stage="pretrain", out_dir=out_dir, data_dir=FIXTURE, seed=seed, device="cpu",
        epochs=2, min_epochs=1, batch_size=32, max_train_rows=120, max_valid_rows=60,
        smoke=True, model_card=False, tahoe_cache=cache, command=["test"],
    )


def _finetune_options(out_dir, seed=5, init="random", init_checkpoint=None,
                      data_dir=None, epochs=2):
    return train_mod.TrainOptions(
        stage="finetune", out_dir=out_dir,
        data_dir=FIXTURE if data_dir is None else data_dir,
        seed=seed, device="cpu", fold=0, epochs=epochs, min_epochs=1, batch_size=32,
        lr=3e-4, init=init, init_checkpoint=init_checkpoint, smoke=True,
        model_card=False, command=["test"],
    )


class EndToEndSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        global FIXTURE, PRETRAIN, RANDOM, WARM, TMP
        TMP = tempfile.mkdtemp(prefix="ptlib-e2e-")
        FIXTURE = _support.fixture_dir()
        root = Path(TMP)
        PRETRAIN = root / "pretrain"
        RANDOM = root / "cm-random"
        WARM = root / "cm-tahoe"
        cls.pretrain_result = train_mod.run_training(_pretrain_options(PRETRAIN))
        cls.random_result = train_mod.run_training(_finetune_options(RANDOM))
        cls.warm_result = train_mod.run_training(
            _finetune_options(
                WARM, init="pretrained",
                init_checkpoint=PRETRAIN / "checkpoint" / "model.pt",
            )
        )

    def test_pretrain_artifacts_and_manifest(self):
        for name in ("manifest.json", "config.json", "fit_history.json",
                     "summary_table.csv", "checkpoint/model.pt"):
            self.assertTrue((PRETRAIN / name).exists(), name)
        manifest = json.loads((PRETRAIN / "manifest.json").read_text())
        self.assertEqual(manifest["kind"], "train")
        self.assertEqual(manifest["stage"], "pretrain")
        self.assertEqual(manifest["status"], "completed")
        self.assertEqual(manifest["structure"]["n_genes"], 18164)
        self.assertEqual(manifest["init"]["mode"], "random")
        self.assertEqual(manifest["metric"]["k"], 200)
        self.assertTrue(manifest["checkpoint"]["reload_identical"])
        self.assertEqual(manifest["noise_policy"]["validation"],
                         "zero noise, model.eval(), fixed row order")
        self.assertIn("sha", manifest["git"])
        self.assertEqual(len(manifest["data"]["artifacts"]), 9)
        self.assertEqual(len(manifest["timings"]["per_epoch_valid_s"]), 2)
        self.assertIsNotNone(manifest["stop_reason"])
        history = json.loads((PRETRAIN / "fit_history.json").read_text())["history"]
        self.assertEqual(len(history), 2)
        self.assertTrue(all(entry["valid_seconds"] > 0 for entry in history))
        self.assertTrue(all(entry["valid_gather_seconds"] >= 0 for entry in history))
        self.assertTrue(all(entry["valid_metrics_seconds"] >= 0 for entry in history))

    def test_checkpoint_payload_and_strict_reload(self):
        payload = checkpoint_mod.read_checkpoint(WARM / "checkpoint" / "model.pt")
        self.assertEqual(set(payload), set(checkpoint_mod.PAYLOAD_KEYS))
        self.assertEqual(payload["n_genes"], 18164)
        model = model_mod.build_model(payload["config"]["structure"])
        checkpoint_mod.load_into_model(model, payload)
        first = model.state_dict()["decoder.recon_decoder.weight"].clone()
        again = checkpoint_mod.read_checkpoint(WARM / "checkpoint" / "model.pt")
        model2 = model_mod.build_model(payload["config"]["structure"])
        checkpoint_mod.load_into_model(model2, again)
        torch.testing.assert_close(
            model2.state_dict()["decoder.recon_decoder.weight"], first
        )

    def test_warm_start_was_strict_and_recorded(self):
        manifest = json.loads((WARM / "manifest.json").read_text())
        self.assertEqual(manifest["init"]["mode"], "pretrained")
        self.assertEqual(
            manifest["init"]["checkpoint_sha256"],
            checkpoint_mod.checkpoint_sha256(PRETRAIN / "checkpoint" / "model.pt"),
        )
        self.assertEqual(manifest["init"]["source_path"],
                         str(PRETRAIN / "checkpoint" / "model.pt"))

    def test_prediction_is_target_free_and_repeatable(self):
        source = data_mod.CMPredictSource(
            data_mod.open_frozen_input(FIXTURE),
            data_mod.FoldTable.load(FIXTURE / "folds.parquet", 5),
            0,
            data_mod.open_frozen_input(FIXTURE).structures(),
        )
        rows = source.rows_for("test")
        payload = checkpoint_mod.read_checkpoint(RANDOM / "checkpoint" / "model.pt")
        struct = checkpoint_mod.structure_from_checkpoint(payload)
        model = model_mod.build_model(struct)
        checkpoint_mod.load_into_model(model, payload)
        model.eval()
        device = torch.device("cpu")
        first = predict_mod.predict_rows(model, source, rows, 10, 32, device)
        second = predict_mod.predict_rows(model, source, rows, 10, 32, device)
        np.testing.assert_array_equal(first[0], second[0])
        np.testing.assert_array_equal(first[1], second[1])
        np.testing.assert_allclose(
            first[1], first[0] - source.baseline(rows), atol=1e-6
        )
        self.assertFalse(hasattr(source, "target"))

    def test_predict_and_evaluate_test_role(self):
        frozen = data_mod.open_frozen_input(FIXTURE)
        pred_dir = Path(TMP) / "pred-test"
        result = predict_mod.run_predict(
            predict_mod.PredictOptions(
                stage="finetune", checkpoint=RANDOM / "checkpoint" / "model.pt",
                out_dir=pred_dir, data_dir=FIXTURE, role="test", fold=0, device="cpu",
                batch_size=32, repeats=2, command=["test"],
            )
        )
        self.assertTrue(result["repeated_inference_identical"])
        manifest = json.loads((pred_dir / "manifest.json").read_text())
        self.assertEqual(manifest["kind"], "predict")
        self.assertGreater(manifest["rows"], 0)
        self.assertTrue(Path(pred_dir / "pred_expr.npy").exists())
        self.assertTrue(Path(pred_dir / "delta.npy").exists())

        eval_dir = Path(TMP) / "eval-test"
        evaluate_mod.run_evaluation(
            evaluate_mod.EvaluateOptions(
                stage="finetune", pred_dir=pred_dir, out_dir=eval_dir,
                data_dir=FIXTURE, role="test", fold=0, command=["test"],
            )
        )
        aggregate = json.loads((eval_dir / "aggregate.json").read_text())
        self.assertEqual(aggregate["headline_k"], 200)
        self.assertGreater(aggregate["n_conditions"], 0)
        macro = aggregate["by_k"]["200"]["macro"]["mean"]
        self.assertTrue(0.0 <= macro <= 1.0)
        self.assertTrue((eval_dir / "per_condition.csv").exists())
        self.assertTrue((eval_dir / "metrics.csv").exists())
        self.assertTrue((eval_dir / "deltas.npz").exists())
        self.assertTrue((eval_dir / "manifest.json").exists())

        card = modelcard.write_model_card(RANDOM, eval_dir=eval_dir)
        text = card.read_text()
        self.assertIn("direction AUPRC@200", text)
        self.assertIn("Deviations (D0–D3)", text)
        self.assertIn("Evaluation", text)

    def test_perfect_prediction_scores_one_and_constant_scores_base_rate(self):
        frozen = data_mod.open_frozen_input(FIXTURE)
        rows, condition_ids, names, _info = evaluate_mod.expected_eval_rows(
            frozen, "finetune", "test", 0
        )
        treated, sorted_rows = h5io.read_layer_rows(frozen.cm_h5ad, rows, frozen.n_genes)
        self.assertTrue(np.array_equal(sorted_rows, rows))
        baseline = data_mod.CMControlPool(frozen).mean_for_rows(rows)
        truth_delta = treated.astype(np.float64) - baseline.astype(np.float64)
        condition = data_mod.CMPredictSource(
            frozen, data_mod.FoldTable.load(FIXTURE / "folds.parquet", 5), 0,
            frozen.structures(),
        )
        meta = condition.meta(rows)
        c_star = np.broadcast_to(condition.c_star, treated.shape)
        # ``pred_expr - c*`` must equal the truth Δ, so a perfect prediction expresses
        # the truth Δ on top of the fold's c* (that is what the model would output).
        perfect_expr = treated.astype(np.float64) - baseline + c_star
        for tag, pred_expr in (
            ("perfect", perfect_expr.astype(np.float32)),
            ("constant", c_star.astype(np.float32)),
        ):
            pred_dir = Path(TMP) / f"pred-{tag}"
            pred_dir.mkdir(parents=True, exist_ok=True)
            if tag == "perfect":
                delta = treated.astype(np.float64) - baseline.astype(np.float64)
            else:
                delta = np.zeros_like(treated, dtype=np.float64)
            np.save(pred_dir / "pred_expr.npy", pred_expr)
            np.save(pred_dir / "delta.npy", delta.astype(np.float32))
            meta.to_csv(pred_dir / "meta.csv", index=False)
            (pred_dir / "manifest.json").write_text(
                json.dumps({"checkpoint": {"sha256": "test"}, "kind": "predict"})
            )
            eval_dir = Path(TMP) / f"eval-{tag}"
            evaluate_mod.run_evaluation(
                evaluate_mod.EvaluateOptions(
                    stage="finetune", pred_dir=pred_dir, out_dir=eval_dir,
                    data_dir=FIXTURE, role="test", fold=0, command=["test"],
                )
            )
            aggregate = json.loads((eval_dir / "aggregate.json").read_text())
            macro = aggregate["by_k"]["200"]["macro"]["mean"]
            if tag == "perfect":
                # float32 storage can break exact ties differently from the float64
                # truth, so a perfect ranking is 1 up to tie-rounding, not bit-exact 1.
                self.assertGreater(macro, 0.999)
            else:
                self.assertLess(macro, 0.2)
                self.assertTrue(np.isfinite(aggregate["regression"]["rmse"]["mean"]))

    def test_pretrain_predict_and_evaluate_valid(self):
        pred_dir = Path(TMP) / "pred-pretrain-valid"
        predict_mod.run_predict(
            predict_mod.PredictOptions(
                stage="pretrain", checkpoint=PRETRAIN / "checkpoint" / "model.pt",
                out_dir=pred_dir, data_dir=FIXTURE, role="valid", device="cpu",
                batch_size=32, repeats=2, command=["test"],
            )
        )
        eval_dir = Path(TMP) / "eval-pretrain-valid"
        evaluate_mod.run_evaluation(
            evaluate_mod.EvaluateOptions(
                stage="pretrain", pred_dir=pred_dir, out_dir=eval_dir,
                data_dir=FIXTURE, role="valid", command=["test"],
            )
        )
        aggregate = json.loads((eval_dir / "aggregate.json").read_text())
        self.assertEqual(aggregate["n_rows"], int(self.pretrain_result.manifest["data_counts"]["n_valid"]))
        macro = aggregate["by_k"]["200"]["macro"]["mean"]
        self.assertTrue(0.0 <= macro <= 1.0)

    def test_extension_points_need_no_pipeline_edits(self):
        """#10 (`--structure`) and #11 (`--train-drugs expanded`) resolve through config only."""
        report = train_mod.check_run(
            train_mod.TrainOptions(
                stage="finetune", out_dir=Path(TMP) / "check-extension", data_dir=FIXTURE,
                device="cpu", fold=0, structure_spec="latent128",
                train_drugs="expanded", model_card=False, command=["test"],
            )
        )
        self.assertEqual(report["model"]["structure"]["z_dim"], 128)
        self.assertEqual(report["model"]["parameter_count"], 7_201_448)
        self.assertGreater(report["counts"]["train_rows_expanded"], 0)
        report256 = train_mod.check_run(
            train_mod.TrainOptions(
                stage="finetune", out_dir=Path(TMP) / "check-extension-256", data_dir=FIXTURE,
                device="cpu", fold=0, structure_spec="hidden256", model_card=False,
                command=["test"],
            )
        )
        self.assertEqual(report256["model"]["parameter_count"], 14_194_280)

    def test_poisoning_test_labels_does_not_change_training(self):
        clean = _support.fixture_dir()
        _, _, test_rows = _support.test_rows_for_fold(clean, 0)
        poisoned = _support.copy_fixture("poisoned-e2e")
        _support.poison_cm_targets(clean, poisoned, test_rows)
        root = Path(TMP)
        clean_run = root / "poison-clean"
        poisoned_run = root / "poison-poisoned"
        train_mod.run_training(
            _finetune_options(clean_run, data_dir=clean, epochs=1)
        )
        train_mod.run_training(
            _finetune_options(poisoned_run, data_dir=poisoned, epochs=1)
        )
        clean_history = json.loads((clean_run / "fit_history.json").read_text())["history"]
        poisoned_history = json.loads(
            (poisoned_run / "fit_history.json").read_text()
        )["history"]
        self.assertEqual(len(clean_history), 1)
        self.assertEqual(
            clean_history[0]["train_loss"], poisoned_history[0]["train_loss"]
        )
        self.assertEqual(
            clean_history[0]["valid_auprc"], poisoned_history[0]["valid_auprc"]
        )
        # ... but the evaluator's truth Δ does change, proving only the evaluator reads it
        clean_pred = root / "poison-clean-pred"
        poisoned_pred = root / "poison-poisoned-pred"
        for run, pred in ((clean_run, clean_pred), (poisoned_run, poisoned_pred)):
            predict_mod.run_predict(
                predict_mod.PredictOptions(
                    stage="finetune", checkpoint=run / "checkpoint" / "model.pt",
                    out_dir=pred, data_dir=clean if "clean" in run.name else poisoned,
                    role="test", fold=0, device="cpu", batch_size=32, repeats=1,
                    command=["test"],
                )
            )
        np.testing.assert_array_equal(
            np.load(clean_pred / "pred_expr.npy"), np.load(poisoned_pred / "pred_expr.npy")
        )
        eval_clean = root / "poison-eval-clean"
        eval_poison = root / "poison-eval-poisoned"
        for pred, out, data_dir in (
            (clean_pred, eval_clean, clean),
            (poisoned_pred, eval_poison, poisoned),
        ):
            evaluate_mod.run_evaluation(
                evaluate_mod.EvaluateOptions(
                    stage="finetune", pred_dir=pred, out_dir=out,
                    data_dir=data_dir, role="test", fold=0, command=["test"],
                )
            )
        a = np.load(eval_clean / "deltas.npz")["true_delta"]
        b = np.load(eval_poison / "deltas.npz")["true_delta"]
        self.assertGreater(float(np.max(np.abs(a - b))), 0.5)


if __name__ == "__main__":
    unittest.main()
