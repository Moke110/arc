"""Locked-config tests: the spec literals, the derived defaults and the digest.

The config is a fixed decision document, not a per-build output: every value is
quoted from the approved spec/ADR or marked `derived_default` with a rationale.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))

import contract


class LockedConfigTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = contract.build_locked_config()

    def test_model_block(self):
        model = self.config["model"]
        self.assertEqual(model["in_dim"], 18164)
        self.assertEqual(model["hidden_layers"], [128])
        self.assertEqual(model["drug_latent"], 64)
        self.assertEqual(model["state_latent"], 64)
        self.assertEqual(model["noise_dim"], 10)
        self.assertEqual(model["adaptor_layers"], [128])
        self.assertEqual(model["drug_dimension"], 1024)
        self.assertFalse(model["signed_delta_head"], "PRnet keeps the expression target")
        self.assertEqual(model["params_expected"], 7_185_000)
        self.assertEqual(model["params_expected"], 386 * 18164 + 173696)

    def test_dose_encoding(self):
        dose = self.config["dose_encoding"]
        self.assertEqual(dose["n_bits"], 1024)
        self.assertEqual(dose["radius"], 2)
        self.assertTrue(dose["use_features"])
        self.assertEqual(dose["dose_transform"], "log10(dose_uM + 1)")
        self.assertEqual(dose["binning"], "none; the two sources' dose grids are not aligned")
        self.assertAlmostEqual(dose["documented_scalings"]["0.2"], 0.07918, places=5)
        self.assertAlmostEqual(dose["documented_scalings"]["10.0"], 1.04139, places=5)

    def test_loss_is_never_clamped(self):
        loss = self.config["loss"]
        self.assertEqual(loss["name"], "Gaussian NLL")
        self.assertFalse(loss["clamp"])
        self.assertIn("never clamp", loss["note"])

    def test_optimisation_recipe(self):
        optim = self.config["optimisation"]
        self.assertEqual(optim["optimizer"], "Adam")
        self.assertEqual(optim["batch_size"], 512)
        self.assertEqual(optim["weight_decay"], 1e-08)
        self.assertEqual(optim["max_epochs"], 100)
        self.assertEqual(optim["min_epochs"], 20)
        self.assertEqual(optim["early_stop_patience"], 20)
        self.assertEqual(optim["plateau_factor"], 0.5)
        self.assertEqual(optim["plateau_patience"], 10)
        self.assertEqual(optim["plateau_mode"], "max",
                         "selection monitors validation AUPRC, which is maximised")
        self.assertEqual(optim["plateau_metric"], "validation direction AUPRC@200")
        self.assertEqual(optim["pretrain_lr"], 0.001)
        self.assertEqual(optim["cm_lr_candidates"], [0.0001, 0.0003, 0.001])
        self.assertIn("full-parameter", optim["finetune"])

    def test_folds_block(self):
        folds = self.config["folds"]
        self.assertEqual(folds["n_folds"], 5)
        self.assertEqual(folds["split_unit"], "drug_key")
        self.assertEqual(folds["seed"], contract.FOLD_SEED)
        self.assertIn("6:2:2", folds["ratios"])
        self.assertIn("NO fold", folds["excluded_drugs"])
        self.assertIn("NOT a plate split", folds["plate_confounding"])
        self.assertEqual(folds["roles"], ["train", "valid", "test"])

    def test_consensus_semantics(self):
        consensus = self.config["consensus_control"]
        self.assertIn("vehicle == 'DMSO'", consensus["candidate_wells"])
        self.assertIn("TRAINING drug", consensus["candidate_wells"])
        self.assertIn("test/valid", consensus["held_out_environment_rule"])
        self.assertIn("equal weight", consensus["cross_plate_weight"])
        self.assertIn("test/valid drug", consensus["forbidden"])
        self.assertIn("fixed per-fold consensus c*",
                      self.config["data_contract"]["cm_validation_test_input"])

    def test_metric_definition(self):
        metric = self.config["metric"]
        self.assertEqual(metric["headline"], "direction AUPRC@200")
        self.assertEqual(metric["column"], "AUPRC@200")
        self.assertEqual(metric["k"], 200)
        self.assertIn("top-k up", metric["true_set"])
        self.assertIn("drug_key x dose_uM", metric["condition"])
        self.assertIn("BEFORE scoring", metric["replicate_rule"])
        self.assertEqual(metric["aggregation"], "macro mean over conditions")
        self.assertTrue(metric["not_a_significance_test"])
        self.assertEqual(metric["selection_metric"], "validation AUPRC@200 (max)")

    def test_validation_block(self):
        validation = self.config["validation"]
        self.assertIn("plate13", validation["tahoe_pretrain"])
        self.assertIn("plate14", validation["tahoe_pretrain"])
        self.assertIn("cell_line", validation["tahoe_matched_control"])
        self.assertIn("c*", validation["cm"])

    def test_repetition_budget(self):
        seeds = self.config["seeds"]
        self.assertEqual(len(seeds["model_seeds"]), 3)
        self.assertTrue(seeds["no_seed_cherry_picking"])
        units = self.config["test_units"]
        self.assertEqual(units["count"], 30)
        self.assertEqual(units["candidates"], 90)
        self.assertEqual(units["factors"], "2 inits x 5 folds x 3 seeds")
        self.assertIn("no seed or fold selection", units["policy"])

    def test_lr_selection_rule(self):
        rule = self.config["lr_selection"]
        self.assertIn("mean over the three seeds", rule["rule"])
        self.assertIn("validation AUPRC@200", rule["rule"])
        self.assertIn("lower LR", rule["tie_break"])

    def test_prediction_policy(self):
        prediction = self.config["prediction"]
        self.assertTrue(prediction["deterministic"])
        self.assertEqual(prediction["validation_test_inference_noise"], "exactly zero")
        self.assertIn("seeded", prediction["train_noise"])

    def test_derived_defaults_are_rationalised(self):
        entries = self.config["derived_defaults"]
        self.assertTrue(entries)
        keys = set()
        for entry in entries:
            self.assertIn("key", entry)
            self.assertIn("value", entry)
            self.assertTrue(entry["rationale"].strip())
            keys.add(entry["key"])
        self.assertEqual(len(keys), len(entries), "derived_default keys are unique")

    def test_digest_is_stable(self):
        first = contract.config_digest(self.config)
        second = contract.config_digest(contract.build_locked_config())
        self.assertEqual(first, second)
        self.assertRegex(first, r"^[0-9a-f]{64}$")

    def test_config_has_no_pending_or_todo(self):
        import json
        text = json.dumps(self.config).lower()
        for marker in ("todo", "tbd", "fixme", "pending"):
            self.assertNotIn(marker, text)


if __name__ == "__main__":
    unittest.main()
