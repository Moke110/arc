"""Locked-config tests: canonical pipeline schema, spec literals, derived defaults.

The top level mirrors the canonical structure/recipe/metric schema ticket #4's
`ptlib.config` reads (structure, recipe, metric, loss, optimizer,
prediction_noise, n_folds, seeds, lrs, dose_encoding) so the pipeline binds to
the frozen values instead of falling back to its own defaults.  The checks are
written independently of #4's code on purpose: a drift between the two must
fail here, not silently default.
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))

import contract

#: ticket #4's `CANONICAL_STRUCTURE_DEFAULTS` / `CANONICAL_RECIPE_DEFAULTS`
PIPELINE_STRUCTURE_KEYS = ("n_genes", "hidden", "z_dim", "c_dim", "n_dim",
                           "drug_dim", "adaptor", "dropout")
PIPELINE_RECIPE_KEYS = ("batch_size", "max_epochs", "min_epochs",
                        "early_stop_patience", "sched_factor", "sched_patience",
                        "min_lr", "weight_decay")


class LockedConfigTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = contract.build_locked_config()

    def test_pipeline_structure_schema(self):
        struct = self.config["structure"]
        for key in PIPELINE_STRUCTURE_KEYS:
            self.assertIn(key, struct)
        self.assertEqual(struct["n_genes"], 18164)
        self.assertEqual(struct["hidden"], [128])
        self.assertEqual(struct["z_dim"], 64)
        self.assertEqual(struct["c_dim"], 64)
        self.assertEqual(struct["n_dim"], 10)
        self.assertEqual(struct["drug_dim"], 1024)
        self.assertEqual(struct["adaptor"], [128])
        self.assertEqual(struct["dropout"], 0.05)

    def test_pipeline_recipe_schema(self):
        recipe = self.config["recipe"]
        for key in PIPELINE_RECIPE_KEYS:
            self.assertIn(key, recipe)
        self.assertEqual(recipe["batch_size"], 512)
        self.assertEqual(recipe["max_epochs"], 100)
        self.assertEqual(recipe["min_epochs"], 20)
        self.assertEqual(recipe["early_stop_patience"], 20)
        self.assertEqual(recipe["sched_factor"], 0.5)
        self.assertEqual(recipe["sched_patience"], 10)
        self.assertEqual(recipe["min_lr"], 1e-08)
        self.assertEqual(recipe["weight_decay"], 1e-08)
        self.assertGreaterEqual(recipe["min_epochs"], 1)
        self.assertLessEqual(recipe["min_epochs"], recipe["max_epochs"])
        self.assertGreater(recipe["batch_size"], 1, "BatchNorm1d needs >1 sample")

    def test_pipeline_metric_schema(self):
        metric = self.config["metric"]
        self.assertEqual(metric["k"], 200)
        self.assertEqual(metric["name"], "direction_auprc")
        self.assertEqual(metric["plateau_mode"], "max")
        self.assertEqual(metric["headline"], "direction AUPRC@200")
        self.assertEqual(metric["column"], "AUPRC@200")
        self.assertIn("drug_key x dose_uM", metric["condition"])
        self.assertIn("BEFORE scoring", metric["replicate_rule"])
        self.assertEqual(metric["aggregation"], "macro mean over conditions")
        self.assertTrue(metric["not_a_significance_test"])
        self.assertEqual(metric["selection_metric"],
                         "validation direction AUPRC@200 (max)")

    def test_pipeline_scalars(self):
        self.assertEqual(self.config["loss"], "gaussian_nll")
        self.assertEqual(self.config["optimizer"], "adam")
        self.assertEqual(self.config["prediction_noise"], "zero")
        self.assertEqual(self.config["n_folds"], 5)
        self.assertEqual(self.config["seeds"], list(contract.MODEL_SEEDS))
        self.assertEqual(len(self.config["seeds"]), 3)
        self.assertEqual(self.config["lrs"]["pretrain"], 0.001)
        self.assertEqual(self.config["lrs"]["cm"], [0.0001, 0.0003, 0.001])

    def test_pipeline_dose_encoding_schema(self):
        dose = self.config["dose_encoding"]
        self.assertEqual(dose["fingerprint"], "FCFP4")
        self.assertEqual(dose["n_bits"], 1024)
        self.assertEqual(dose["radius"], 2)
        self.assertTrue(dose["use_features"])
        self.assertEqual(dose["dose_scale"], "log10(dose_uM+1)")
        self.assertEqual(dose["binning"], "none; the two sources' dose grids are not aligned")
        self.assertAlmostEqual(dose["documented_scalings"]["0.2"], 0.07918, places=5)
        self.assertAlmostEqual(dose["documented_scalings"]["10.0"], 1.04139, places=5)

    def test_gaussian_nll_never_clamped(self):
        detail = self.config["loss_detail"]
        self.assertFalse(detail["clamp"])
        self.assertIn("never clamp", detail["note"])
        self.assertIn("gaussian_nll_loss", detail["implementation"])

    def test_model_detail(self):
        detail = self.config["model_detail"]
        self.assertEqual(detail["params_expected"], 7_185_000)
        self.assertEqual(detail["params_expected"], 386 * 18164 + 173696)
        self.assertEqual(detail["in_dim"], 18164)
        self.assertEqual(detail["comb_num"], 1)
        self.assertFalse(detail["signed_delta_head"])
        self.assertIn("full-parameter", detail["finetune"])
        self.assertEqual(detail["plateau_metric"], "validation direction AUPRC@200")
        self.assertEqual(detail["plateau_mode"], "max")

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

    def test_validation_block(self):
        validation = self.config["validation"]
        self.assertIn("plate13", validation["tahoe_pretrain"])
        self.assertIn("plate14", validation["tahoe_pretrain"])
        self.assertIn("cell_line", validation["tahoe_matched_control"])
        self.assertIn("c*", validation["cm"])

    def test_repetition_budget(self):
        seeds = self.config["seeds_detail"]
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
        text = json.dumps(self.config).lower()
        for marker in ("todo", "tbd", "fixme", "pending"):
            self.assertNotIn(marker, text)


if __name__ == "__main__":
    unittest.main()
